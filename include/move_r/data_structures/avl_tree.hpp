#pragma once

#include <functional>

/**
 * @brief balanced binary search tree, for each node |height of left child - height of right child| <= 1 holds
 * @tparam T value type
 * @tparam Compare comparison class
 */
template <typename T, typename Compare = std::less<T>>
class avl_tree {
    public:

    /**
     * @brief node in an avl_tree
     * @tparam T value type
     */
    struct avl_node {
        T v; // value
        avl_node *p; // parent
        avl_node *lc; // left child
        avl_node *rc; // right child
        uint8_t h; // height

        /**
         * @brief creates an avl_node with parent p
         */
        avl_node(avl_node* p) {
            this->p = p;
        }

        /**
         * @brief creates a copy of another avl_node
         * @param other the avl_node to copy
         */
        void copy_from_other(const avl_node& other) {
            this->v = other.v;
            if (other.lc != NULL) {
                this->lc = new avl_node(this);
                this->lc->copy_from_other(*(other.lc));
                this->lc->p = this;
            } else {
                this->lc = NULL;
            }
            if (other.rc != NULL) {
                this->rc = new avl_node(this);
                this->rc->copy_from_other(*(other.rc));
                this->rc->p = this;
            } else {
                this->rc = NULL;
            }
            this->h = std::max(ht(this->lc),ht(this->rc))+1;
        }

        /**
         * @brief moves an avl_node into this avl_node
         * @param other the avl_node to move
         */
        void move_from_other(avl_node&& other) {
            this->v = std::move(other.v);
            this->p = other.p;
            this->lc = other.lc;
            this->rc = other.rc;
            this->h = other.h;

            other.v = T();
            other.p = NULL;
            other.lc = NULL;
            other.rc = NULL;
            other.h = 0;

            if (lc != NULL) lc->p = this;
            if (rc != NULL) rc->p = this;

            if (p != NULL)  {
                if (&other == p->lc) {
                    p->lc = this;
                } else {
                    p->rc = this;
                }
            }
        }

        /**
         * @brief returns the hight of the node n if n != NULL, else returns 0
         * @param n an avl_node in the avl_tree
         * @return height of n
         */
        inline static uint8_t ht(avl_node *n) {
            return n != NULL ? n->h : 0;
        }
        
        /**
         * @brief recursively deletes all nodes in the subtree of node n
         * @param n an avl_node
         */
        inline void delete_subtree() {
            this->p = NULL;
            delete this;
        }

        public:
        avl_node() = default;
        avl_node(avl_node&& other) {move_from_other(std::move(other));}
        avl_node(const avl_node& other) {copy_from_other(other);}
        avl_node& operator=(avl_node&& other) {move_from_other(std::move(other));return *this;}
        avl_node& operator=(const avl_node& other) {copy_from_other(other);return *this;}
        ~avl_node() {}

        /**
         * @brief creates an avl_node with value v
         */
        avl_node(T v) {
            this->v = v;
            lc = rc = p = NULL;
            h = 0;
        }

        /**
         * @brief returns the next avl_node in the avl_tree
         * @return the next avl_node in the avl_tree
         */
        inline avl_node* nxt() {
            avl_node *cur = this;
            if (cur->rc != NULL) {
                cur = cur->rc;
                while (cur->lc != NULL) {
                    cur = cur->lc;
                }
            } else {
                while (cur->p != NULL && cur == cur->p->rc) {
                    cur = cur->p;
                }
                cur = cur->p;
            }
            return cur;
        }

        /**
         * @brief returns the previous avl_node in the avl_tree
         * @return the previous avl_node in the avl_tree
         */
        inline avl_node* prv() {
            avl_node *cur = this;
            if (cur->lc != NULL) {
                cur = cur->lc;
                while (cur->rc != NULL) {
                    cur = cur->rc;
                }
            } else {
                while (cur->p != NULL && cur == cur->p->lc) {
                    cur = cur->p;
                }
                cur = cur->p;
            }
            return cur;
        }

        /**
         * @brief checks if the avl_node is a leaf
         * @return whether the avl_node is a leaf
         */
        bool is_leaf() {
            return lc == NULL && rc == NULL;
        }
    };

    protected:
    avl_node *r = NULL; // root of the avl_tree
    avl_node *fst = NULL; // first node (node with the smallest value)
    avl_node *lst = NULL; // last node (node with the greatest value)
    uint64_t s = 0; // size
    uint8_t h = 0; // height

    static bool eq(const T& v1, const T& v2) {return !Compare()(v1,v2) && !Compare()(v2,v1);}; // comparison function "equals" on values of type T
    static bool lt(const T& v1, const T& v2) {return Compare()(v1,v2);}; // comparison function "less than" on values of type T
    static bool gt(const T& v1, const T& v2) {return Compare()(v2,v1);}; // comparison function "greater than" on values of type T
    static bool leq(const T& v1, const T& v2) {return Compare()(v1,v2) || !Compare()(v2,v1);}; // comparison function "less than or equal to" on values of type T
    static bool geq(const T& v1, const T& v2) {return Compare()(v2,v1) || !Compare()(v1,v2);}; // comparison function "greater than or equal to" on values of type T

    /**
     * @brief returns the hight of the node n if n != NULL, else returns 0
     * @param n an avl_node in the avl_tree
     * @return height of n
     */
    inline static uint8_t ht(avl_node *n) {
        return avl_node::ht(n);
    }

    /**
     * @brief updates the height of the subtree of node n
     * @param n an avl_node in the avl_tree
     * @return whether the height of n changed
     */
    inline static bool update_height(avl_node *n) {
        uint8_t n_h = n->h;
        n->h = std::max(ht(n->lc),ht(n->rc))+1;
        return n->h != n_h;
    }

    /**
     * @brief left rotation around node x
     * @param x an avl_node in the avl_tree, must have a right child
     */
    inline void rotate_left(avl_node *x) {
        avl_node *y = x->rc;
        x->rc = y->lc;
        if (y->lc != NULL) {
            y->lc->p = x;
        }
        y->p = x->p;
        if (x->p == NULL) {
            r = y;
        } else if (x == x->p->lc) {
            x->p->lc = y;
        } else {
            x->p->rc = y;
        }
        y->lc = x;
        x->p = y;
        update_height(x);
        update_height(y);
    }

    /**
     * @brief right rotation around node y
     * @param y an avl_node in the avl_tree, must have a left child
     */
    inline void rotate_right(avl_node *y) {
        avl_node *x = y->lc;
        y->lc = x->rc;
        if (x->rc != NULL) {
            x->rc->p = y;
        }
        x->p = y->p;
        if (y->p == NULL) {
            r = x;
        } else if (y == y->p->lc) {
            y->p->lc = x;
        } else {
            y->p->rc = x;
        }
        x->rc = y;
        y->p = x;
        update_height(y);
        update_height(x);
    }

    /**
     * @brief balances the node n
     * @param n an avl_node in the avl_tree
     * @return whether n was a-heavy
     */
    inline bool balance(avl_node *n) {
        if (ht(n->lc) > ht(n->rc)+1) {
            if (ht(n->lc->lc) < ht(n->lc->rc)) {
                rotate_left(n->lc);
            }
            rotate_right(n);
        } else if (ht(n->rc) > ht(n->lc)+1) {
            if (ht(n->rc->rc) < ht(n->rc->lc)) {
                rotate_right(n->rc);
            }
            rotate_left(n);
        } else {
            return false;
        }
        return true;
    }

    /**
     * @brief balances all nodes starting from nf up to nt
     * @param nf an avl_node in the avl_tree, has to be in the subtree of nt
     * @param nt an avl_node in the avl_tree
     */
    inline void balance_from_to(avl_node *nf, avl_node *nt) {
        while (nf != nt) {
            if (nf == nf->p->rc) {
                nf = nf->p;
                if (!update_height(nf->rc) & !balance(nf->rc)) return;
            } else {
                nf = nf->p;
                if (!update_height(nf->lc) & !balance(nf->lc)) return;
            }
        }
        update_height(nt);
        balance(nt);
    }

    /**
     * @brief returns the node with the smallest value in the subtree of node n
     * @param n an avl_node in the avl_tree
     * @return the node with the smallest value in the subtree of node n
     */
    inline static avl_node* min(avl_node *n) {
        while (n->lc != NULL) {
            n = n->lc;
        }
        return n;
    }

    /**
     * @brief returns the node with the greatest value in the subtree of node n
     * @param n an avl_node in the avl_tree
     * @return the node with the greatest value in the subtree of node n
     */
    inline static avl_node* max(avl_node *n) {
        while (n->rc != NULL) {
            n = n->rc;
        }
        return n;
    }

    /**
     * @brief removes the node nr in the subtree of node n
     * @param n_rem an avl_node in the avl_tree, must be in the subtree of node n
     * @param n an avl_node in the avl_tree
     */
    inline void remove_node_in(avl_node *n_rem, avl_node *n) {
        avl_node *n_rem_p = n_rem->p;
        if (n_rem->lc == NULL) {
            s--;
            if (n_rem->rc == NULL) {
                if (n_rem->p != NULL) {
                    if (n_rem == n_rem->p->lc) {
                        n_rem->p->lc = NULL;
                    } else {
                        n_rem->p->rc = NULL;
                    }
                } else {
                    r = NULL;
                }
            } else {
                if (n_rem->p != NULL) {
                    if (n_rem == n_rem->p->lc) {
                        n_rem->p->lc = n_rem->rc;
                    } else {
                        n_rem->p->rc = n_rem->rc;
                    }
                    n_rem->rc->p = n_rem->p;
                } else {
                    n_rem->rc->p = NULL;
                    r = n_rem->rc;
                }
                n_rem->rc = NULL;
            }
        } else if (n_rem->rc == NULL) {
            s--;
            if (n_rem->p != NULL) {
                if (n_rem == n_rem->p->lc) {
                    n_rem->p->lc = n_rem->lc;
                } else {
                    n_rem->p->rc = n_rem->lc;
                }
                n_rem->lc->p = n_rem->p;
            } else {
                n_rem->lc->p = NULL;
                r = n_rem->lc;
            }
            n_rem->lc = NULL;
        } else {
            avl_node *n_min = min(n_rem->rc);
            avl_node n_tmp = *n_min;

            n_min->p = n_rem->p;
            if (n_rem->p != NULL) {
                if (n_rem == n_rem->p->lc) {
                    n_rem->p->lc = n_min;
                } else {
                    n_rem->p->rc = n_min;
                }
            }
            n_min->lc = n_rem->lc;
            n_rem->lc->p = n_min;
            n_min->h = n_rem->h;
            if (n_min->rc != NULL) {
                n_min->rc->p = n_rem;
            }
            if (n_min == n_rem->rc) {
                n_min->rc = n_rem;
                n_rem->p = n_min;
            } else {
                n_min->rc = n_rem->rc;
                n_rem->rc->p = n_min;
                n_rem->p = n_tmp.p;
                n_tmp.p->lc = n_rem;
            }
            n_rem->rc = n_tmp.rc;
            n_rem->lc = NULL;
            n_rem->h = n_tmp.h;

            remove_node_in(n_rem,n_min->rc);
        }
        if (n_rem_p != NULL) {
            balance_from_to(n_rem_p,n);
        }
    }

    /**
     * @brief creates a balanced avl subtree of the nodes in at(l),at(l+1),...,at(r) and returns it's root node
     * @param l l in [0..|nds|-1]
     * @param r r in [0..|nds|-1], l <= r
     * @param at function returning the node at a given position
     * @param max_tasks max number of tasks to start
     * @return the root of the avl subtree
     */
    template <typename uint_t>
    inline static avl_node* build_subtree(uint_t l, uint_t r, std::function<avl_node*(uint_t)> &at, uint16_t max_tasks = 1) {
        if (r == l) {
            avl_node* n_l = at(l);
            n_l->rc = n_l->lc = NULL;
            n_l->h = 0;
            return n_l;
        } else if (r == l+1) {
            avl_node* n_l = at(l);
            avl_node* n_r = at(r);
            n_r->rc = n_l->lc = n_l->rc = NULL;
            n_r->lc = n_l;
            n_r->h = 1;
            n_l->p = n_r;
            n_l->h = 0;
            return n_r;
        } else {
            uint_t m = l+(r-l)/2;
            avl_node* n_m = at(m);
            if (max_tasks > 1) {
                #pragma omp task
                {
                    n_m->lc = build_subtree(l,m-1,at,max_tasks/2);
                }
                n_m->rc = build_subtree(m+1,r,at,max_tasks/2);
                #pragma omp taskwait
            } else {
                n_m->lc = build_subtree(l,m-1,at,1);
                n_m->rc = build_subtree(m+1,r,at,1);
            }
            n_m->lc->p = n_m;
            n_m->rc->p = n_m;
            n_m->h = std::max(n_m->lc->h,n_m->rc->h)+1;
            return n_m;
        }
    }

    /**
     * @brief creates a copy of another avl_tree
     * @param other the avl_tree to copy
     */
    void copy_from_other(const avl_tree& other) {
        if (other.s > 0) {
            this->s = other.s;
            this->h = other.h;
            this->r = new avl_node();
            *(this->r) = *(other.r);
            this->r->p = NULL;
            this->fst = min(r);
            this->lst = max(r);
        }
    }

    /**
     * @brief moves an avl_tree into this avl_tree
     * @param other the avl_tree to move
     */
    void move_from_other(avl_tree&& other) {
        this->r = other.r;
        this->fst = other.fst;
        this->lst = other.lst;
        this->s = other.s;
        this->h = other.h;

        other.r = NULL;
        other.fst = NULL;
        other.lst = NULL;
        other.s = 0;
        other.h = 0;
    }

    public:
    avl_tree() = default;
    avl_tree(avl_tree&& other) {move_from_other(std::move(other));}
    avl_tree(const avl_tree& other) {copy_from_other(other);}
    avl_tree& operator=(avl_tree&& other) {move_from_other(std::move(other));return *this;}
    avl_tree& operator=(const avl_tree& other) {copy_from_other(other);return *this;}

    /**
     * @brief deletes the avl_tree but not it's nodes
     */
    ~avl_tree() {
        if (r != NULL) r->delete_subtree();
        fst = lst = r = NULL;
    }

    /**
     * @brief returns the height of the avl_tree
     * @return height of the avl_tree
     */
    inline uint8_t height() {
        return h;
    }

    /**
     * @brief returns the number of elements in the avl_tree
     * @return number of elements in the avl_tree
     */
    inline uint64_t size() {
        return s;
    }

    /**
     * @brief returns whether the avl_tree is empty
     * @return whether the avl_tree is empty
     */
    inline bool empty() {
        return s == 0;
    }

    /**
     * @brief deletes all nodes in the avl_tree
     */
    inline void delete_nodes(){
        if (!empty()) {
            r->delete_subtree();
            fst = lst = r = NULL;
            h = s = 0;
        }
    }

    /**
     * @brief disconnects all nodes in the avl_tree
     */
    inline void disconnect_nodes() {
        fst = lst = r = NULL;
        h = s = 0;
    }

    /**
     * @brief inserts the nodes in at(l),at(l+1),...,at(r) into the avl_tree if it is empty
     * @param l l in [0..|nds|-1]
     * @param r r in [0..|nds|-1], l <= r
     * @param at function returning the node at a given position
     * @param max_tasks max number of tasks to start
     */
    template <typename uint_t>
    void insert_array(uint_t l, uint_t r, std::function<avl_node*(uint_t)> &at, uint16_t max_tasks = 1) {
        if (empty() && l >= 0 && r >= l) {
            this->r = build_subtree(l,r,at,omp_in_parallel() ? max_tasks : 1);
            this->r->p = NULL;
            this->fst = at(l);
            this->lst = at(r);
            h = this->r->h;
            s = r-l+1;
        }
    }

    /**
     * @brief returns the node with the smallest value in the avl_tree
     * @return the node with the smallest value in the avl_tree if the avl_tree is not empty, else NULL
     */
    inline avl_node* min() {
        return fst;
    }

    /**
     * @brief returns the node with the second smallest value in the avl_tree
     * @return the node with the second smallest value in the avl_tree if the avl_tree has at least 2 nodes, else NULL
     */
    inline avl_node* second_smallest() {
        if (fst->rc != NULL) {
            if (fst->rc->lc != NULL) {
                return fst->rc->lc;
            } else {
                return fst->rc;
            }
        } else {
            return fst->p;
        }
    }

    /**
     * @brief returns the node with the greatest value in the avl_tree
     * @return the node with the greatest value in the avl_tree if the avl_tree is not empty, else NULL
     */
    inline avl_node* max() {
        return lst;
    }

    /**
     * @brief returns the node with the second largest value in the avl_tree
     * @return the node with the second largest value in the avl_tree if the avl_tree has at least 2 nodes, else NULL
     */
    inline avl_node* second_largest() {
        if (lst->lc != NULL) {
            if (lst->lc->rc != NULL) {
                return lst->lc->rc;
            } else {
                return lst->lc;
            }
        } else {
            return lst->p;
        }
    }

    /**
     * @brief creates and inserts a node with the value v into the avl_tree or
     *        updates the node in the avl_tree with a value equal to v
     * @param v value
     * @param n an avl_node in the avl_tree, v must be able to be inserted into the subtree of node n
     * @return the node in the avl_tree with the value v
     */
    inline avl_node* emplace_hint(T &&v, avl_node *n) {
        if (empty()) {
            r = new avl_node(v);
            h = s = 1;
            fst = lst = r;
            return r;
        } else {
            n = find(v,n);
            if (eq(n->v,v)) {
                n->v = v;
                return n;
            } else {
                avl_node *n_new = new avl_node(v);
                if (lt(n_new->v,n->v)) {
                    n->lc = n_new;
                    n_new->p = n;
                } else {
                    n->rc = n_new;
                    n_new->p = n;
                }
                balance_from_to(n_new->p,r);
                h = r->h;
                s++;
                if (lt(v,fst->v)) {
                    fst = n_new;
                } else if (gt(v,lst->v)) {
                    lst = n_new;
                }
                return n_new;
            }
        }
    }

    /**
     * @brief creates and inserts a node with the value v into the avl_tree or
     *        updates the node in the avl_tree with a value equal to v
     * @param v value
     * @param n an avl_node in the avl_tree, v must be able to be inserted into the subtree of node n
     * @return the node in the avl_tree with the value v
     */
    inline avl_node* emplace_hint(T &v, avl_node *n) {
        return emplace_hint(std::move(v),n);
    }

    /**
     * @brief creates and inserts a node with the value v into the avl_tree or
     *        updates the node in the avl_tree with a value equal to v
     * @param v value
     * @return the node in the avl_tree with the value v
     */
    inline avl_node* insert(T &&v) {
        return emplace_hint(v,r);
    }

    /**
     * @brief creates and inserts a node with the value v into the avl_tree or
     *        updates the node in the avl_tree with a value equal to v
     * @param v value
     * @return the node in the avl_tree with the value v
     */
    inline avl_node* insert(T &v) {
        return emplace_hint(std::move(v),r);
    }

    /**
     * @brief inserts the node n into the avl_tree
     * @param n an avl_node, it and it's value must not be in the avl_tree
     * @param nin an avl_node in the avl_tree, n must be able to be inserted into the subtree of node na
     * @return n
     */
    inline avl_node* insert_hint(avl_node *n, avl_node *n_in) {
        avl_node *n_at = find(n->v,n_in);
        if (lt(n->v,n_at->v)) {
            n_at->lc = n;
            n->p = n_at;
        } else {
            n_at->rc = n;
            n->p = n_at;
        }
        balance_from_to(n->p,r);
        h = r->h;
        s++;
        if (lt(n->v,fst->v)) {
            fst = n;
        } else if (gt(n->v,lst->v)) {
            lst = n;
        }
        return n;
    }

    /**
     * @brief inserts the node n into the avl_tree
     * @param n an avl_node, it and it's value must not be in the avl_tree
     * @return n
     */
    inline avl_node* insert_node(avl_node *n) {
        if (empty()) {
            fst = lst = r = n;
            h = s = 1;
            return r;
        } else {
            return insert_hint(n,r);
        }
    }

    /**
     * @brief removes node n from the avl_tree
     * @param n an avl_node in the avl_tree
     * @return the node that has been removed
     */
    inline avl_node* remove_node(avl_node *n) {
        if (s == 1) {
            r = fst = lst = NULL;
            h = s = 0;
        } else {
            if (n == fst) {
                fst = second_smallest();
            } else if (n == lst) {
                lst = second_largest();
            }
            remove_node_in(n,r);
            h = r->h;
        }
        return n;
    }

    /**
     * @brief removes the node with a value equal to v from the avl_tree
     * @param v value
     * @return the node that has been removed if it was in the avl_tree
     */
    inline avl_node* remove(T &&v) {
        avl_node *n = find(v);
        if (n == NULL || !eq(n->v,v)) return NULL;
        return remove_node(n);
    }

    /**
     * @brief removes the node with a value equal to v from the avl_tree
     * @param v value
     * @return whether there was a node with a value equal to v in the avl_tree
     */
    inline avl_node* remove(T &v) {
        return remove(std::move(v));
    }

    /**
     * @brief searches for a node with value v in the avl_tree until it was found or a leaf has been reached
     * @param v value
     * @param n an avl_node in the avl_tree, at which the search starts
     * @return the node with a value equal to v or a leaf at which a node with value v can be inserted
     *         if the avl_tree is not empty, else NULL
     */
    inline avl_node* find(T &&v, avl_node *n) {
        if (empty()) return NULL;
        while (true) {
            if (gt(v,n->v)) {
                if (n->rc != NULL) {
                    n = n->rc;
                } else {
                    break;
                }
            } else if (lt(v,n->v)) {
                if (n->lc != NULL) {
                    n = n->lc;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        return n;
    }

    /**
     * @brief searches for a node with value v in the avl_tree until it was found or a leaf has been reached
     * @param v value
     * @param n an avl_node in the avl_tree, at which the search starts
     * @return the node with a value equal to v or a leaf at which a node with value v can be inserted
     *         if the avl_tree is not empty, else NULL
     */
    inline avl_node* find(T &v, avl_node *n) {
        return find(std::move(v),n);
    }

    /**
     * @brief searches for a node with value v in the avl_tree until it was found or a leaf has been reached
     * @param v value
     * @return avl_node* the node with a value equal to v or a leaf at which a node with value v can be inserted
     *         if the avl_tree is not empty, else NULL
     */
    inline avl_node* find(T &&v) {
        return find(v,r);
    }

    /**
     * @brief searches for a node with value v in the avl_tree until it was found or a leaf has been reached
     * @param v value
     * @return avl_node* the node with a value equal to v or a leaf at which a node with value v can be inserted
     *         if the avl_tree is not empty, else NULL
     */
    inline avl_node* find(T &v) {
        return find(std::move(v),r);
    }

    /**
     * @brief returns the node with the smallest value greater than or equal to v
     * @param v value
     * @return avl_node* the node with the smallest value greater than or equal to v, if it exists
     *         if not all nodes' values are smaller than v, else NULL 
     */
    inline avl_node* min_geq(T &&v) {
        if (empty()) return NULL;
        avl_node *n = r;
        avl_node *min = NULL;
        while (true) {
            if (lt(n->v,v)) {
                if (n->rc != NULL) {
                    n = n->rc;
                } else {
                    break;
                }
            } else if (eq(n->v,v)) {
                return n;
            } else if (n->lc != NULL) {
                if (lt(n->lc->v,v)) {
                    min = n;
                }
                n = n->lc;
            } else {
                break;
            }
        }
        if (geq(n->v,v)) {
            return n;
        }
        return min;
    }

    /**
     * @brief returns the node with the smallest value greater than or equal to v
     * @param v value
     * @return avl_node* the node with the smallest value greater than or equal to v, if it exists
     *         if not all nodes' values are smaller than v, else NULL 
     */
    inline avl_node* min_geq(T &v) {
        return min_geq(std::move(v));
    };

    /**
     * @brief returns the node with the greatest value less than or equal to v
     * @param v value
     * @return avl_node* the node with the greatest value less than or equal to v, if it exists
     *         if not all nodes' values are greater than v, else NULL 
     */
    inline avl_node* max_leq(T &&v) {
        if (empty()) return NULL;
        avl_node *n = r;
        avl_node* max = NULL;
        while (true) {
            if (gt(n->v,v)) {
                if (n->lc != NULL) {
                    n = n->lc;
                } else {
                    break;
                }
            } else if (eq(n->v,v)) {
                return n;
            } else if (n->rc != NULL) {
                if (gt(n->rc->v,v)) {
                    max = n;
                }
                n = n->rc;
            } else {
                break;
            }
        }
        if (leq(n->v,v)) {
            return n;
        }
        return max;
    }

    /**
     * @brief returns the node with the greatest value less than or equal to v
     * @param v value
     * @return avl_node* the node with the greatest value less than or equal to v, if it exists
     *         if not all nodes' values are greater than v, else NULL 
     */
    inline avl_node* max_leq(T &v) {
        return max_leq(std::move(v));
    }

    /**
     * @brief iterator for an avl_tree
     */
    class avl_it {
        protected:
        avl_tree *t; // the avl_tree, the iterator iterates through
        avl_node *cur; // the node the iterator points to

        public:
        /**
         * @brief creates an avl_it pointing to the node n in the avl_tree t
         * @param t an avl_tree
         * @param n an avl_node in t
         */
        avl_it(avl_tree *t, avl_node *n) {
            this->t = t;
            this->cur = n;
        }

        /**
         * @brief deletes the avl_it
         */
        ~avl_it() {
            t = NULL;
            cur = NULL;
        }
        
        /**
         * @brief checks if the iterator can iterate forward
         * @return whether it can iterate forward
         */
        inline bool has_next() {
            return avl_tree::lt(cur->v,t->lst->v);
        }

        /**
         * @brief checks if the iterator can iterate backward
         * @return whether it can iterate backward
         */
        inline bool has_prev() {
            return avl_tree::gt(cur->v,t->fst->v);
        }

        /**
         * @brief returns the value of the node, the iterator points to
         * @return the node, the iterator points to
         */
        inline avl_node* current() {
            return cur;
        }

        /**
         * @brief iterates forward, has_next() must return true
         * @return the node the iterator points to after iterating forward
         */
        inline avl_node* next() {
            cur = cur->nxt();
            return cur;
        }

        /**
         * @brief iterates forward, has_pred() must return true
         * @return the node the iterator points to after iterating backward
         */
        inline avl_node* previous() {
            cur = cur->prv();
            return cur;
        }

        /**
         * @brief points the iterator to the node n
         * @param n an avl_node in t
         */
        inline void set(avl_node *n) {
            cur = n;
        }
    };

    /**
     * @brief returns an iterator pointing to the node n
     * @param n an avl_node in the avl_tree
     * @return an iterator
     */
    inline avl_tree::avl_it iterator(avl_node *n) {
        return avl_tree::avl_it(this,n);
    }

    /**
     * @brief returns an iterator pointing to the min of the avl_tree if it is not empty
     * @return an iterator
     */
    inline avl_tree::avl_it iterator() {
        return avl_tree::avl_it(this,fst);
    }
};