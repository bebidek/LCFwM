#ifndef SUFFIX_TREE_HPP
#define SUFFIX_TREE_HPP

#include <string_view>
#include <vector>
#include <functional>

struct suffix_tree {
    // fields
    struct node_t {
        node_t *parent, *child, *sibling, *link;
        int start, depth;
    };

    std::string_view text;
    node_t *root;
    int size;


    // helper functions
    node_t** find_child_ptr(node_t *node, char x) {
        for (auto child_ptr = &node->child; *child_ptr; child_ptr = &(*child_ptr)->sibling)
            if (text[(*child_ptr)->start + node->depth] == x)
                return child_ptr;
        return nullptr;
    }

    node_t* split_edge(node_t **edge, int depth) {
        node_t *down_node = *edge;
        node_t *up_node = down_node->parent;

        node_t *mid_node = new node_t;
        mid_node->parent = up_node;
        mid_node->child = down_node;
        mid_node->sibling = down_node->sibling;
        mid_node->link = nullptr;
        mid_node->start = down_node->start;
        mid_node->depth = depth;

        *edge = mid_node;
        down_node->parent = mid_node;
        down_node->sibling = nullptr;

        size++;
        return mid_node;
    }

    node_t* create_leaf(node_t *parent, int start) {
        node_t *leaf = new node_t;
        
        leaf->parent = parent;
        leaf->child = nullptr;
        leaf->sibling = parent->child;
        leaf->link = nullptr;
        leaf->start = start;
        leaf->depth = text.length() - start;

        parent->child = leaf;

        size++;
        return leaf;
    }


    // construction / destruction algorithm
    suffix_tree(std::string_view text) : text(text) {
        // create root vertex
        root = new node_t;
        root->parent = root->child = root->sibling = nullptr;
        root->link = root;
        root->start = root->depth = 0;
        size = 1;

        // iterate over suffixes
        node_t *current_node = root;
        int current_depth = 0;

        for (int i = 0; i < (int)text.length(); i++) {
            // go as far as you can, add leaf vertex
            node_t **next_node_ptr;
            while ((next_node_ptr = find_child_ptr(current_node, text[i + current_depth]))) {
                current_depth++;
                while (current_depth < (*next_node_ptr)->depth)
                    if (text[i + current_depth] == text[(*next_node_ptr)->start + current_depth])
                        current_depth++;
                    else {
                        split_edge(next_node_ptr, current_depth);
                        break;
                    }
                current_node = *next_node_ptr;
            }
            create_leaf(current_node, i);

            // update links
            if (!current_node->link) {
                node_t **link_node_ptr = &current_node->parent->link; // will be overwritten at least once, so top-level pointer isn't relevant
                while ((*link_node_ptr)->depth < current_depth - 1)
                    link_node_ptr = find_child_ptr(*link_node_ptr, text[i + 1 + (*link_node_ptr)->depth]);
                if ((*link_node_ptr)->depth > current_depth - 1)
                    current_node->link = split_edge(link_node_ptr, current_depth - 1);
                else
                    current_node->link = *link_node_ptr;
            }

            current_node = current_node->link;
            current_depth = current_node->depth;
        }
    }

    suffix_tree(const suffix_tree&) = delete;

    ~suffix_tree() {
        std::function<void(node_t*)> del = [&del](node_t *node) {
            if (node->child)
                del(node->child);
            if (node->sibling)
                del(node->sibling);
            delete node;
        };
        del(root);
    }
};

struct st_lca {
    typedef suffix_tree::node_t node_t;

    // fields
    std::vector<int> leaf_pos; // start -> position in tab[0]
    std::vector<node_t**> tabs;

    // helpers
    static node_t* min_node(node_t *node1, node_t *node2) {
        return node1->depth < node2->depth ? node1 : node2;
    }

    // construction algorithm
    st_lca(suffix_tree& tree) {
        int size = 2 * tree.size - 1;
        
        // first row - euler cycle
        tabs.push_back(new node_t*[size]);
        leaf_pos.resize(tree.text.length());
        int it = 0;
        
        std::function<void(node_t*)> euler = [&](node_t* node) {
            if (!node->child)
                leaf_pos[node->start] = it;
            tabs[0][it++] = node;

            for (auto child = node->child; child; child = child->sibling) {
                euler(child);
                tabs[0][it++] = node;
            }
        };
        euler(tree.root);

        // next rows - combining information from the previous ones
        for (int s = 2; s < size; s *= 2) {
            node_t** current_row = new node_t*[size];
            for (int i = 0; i < size - s + 1; i++)
                current_row[i] = min_node(tabs.back()[i], tabs.back()[i + s / 2]);
            tabs.push_back(current_row);
        }
    }

    ~st_lca() {
        for (auto row : tabs)
            delete []row;
    }

    // queries
    node_t* lca(int start1, int start2) {
        int pos1 = leaf_pos[start1];
        int pos2 = leaf_pos[start2];

        if (pos1 == pos2)
            return tabs[0][pos1];
        if (pos1 > pos2)
            std::swap(pos1, pos2);

        int row = 64 - __builtin_clzll(pos2 - pos1) - 1;
        return min_node(tabs[row][pos1], tabs[row][pos2 - (1<<row) + 1]);
    }
};

#endif
