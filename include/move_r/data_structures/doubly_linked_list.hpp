#pragma once

/**
 * @brief doubly linked list
 * @tparam T value type
 */
template <typename T>
class doubly_linked_list {
    public:
    /**
     * @brief node in a dll_list
     */
    struct doubly_linked_list_node {
        T v; // value
        doubly_linked_list_node *pr; // predecesor
        doubly_linked_list_node *sc; // successor

        /**
         * @brief creates a copy of another doubly_linked_list_node
         * @param other the doubly_linked_list_node to copy
         */
        void copy_from_other(const doubly_linked_list_node& other) {
            this->v = other.v;
            this->pr = NULL;
            this->sc = NULL;
        }

        /**
         * @brief moves an doubly_linked_list_node into this doubly_linked_list_node
         * @param other the doubly_linked_list_node to move
         */
        void move_from_other(doubly_linked_list_node&& other) {
            this->v = std::move(other.v);
            this->pr = other.pr;
            this->sc = other.sc;

            other.v = T();
            other.pr = NULL;
            other.sc = NULL;

            if (pr != NULL) pr->sc = this;
            if (sc != NULL) sc->pr = this;
        }
        
        public:
        doubly_linked_list_node() = default;
        doubly_linked_list_node(doubly_linked_list_node&& other) {move_from_other(std::move(other));}
        doubly_linked_list_node(const doubly_linked_list_node& other) {copy_from_other(other);}
        doubly_linked_list_node& operator=(doubly_linked_list_node&& other) {move_from_other(std::move(other));return *this;}
        doubly_linked_list_node& operator=(const doubly_linked_list_node& other) {copy_from_other(other);return *this;}
        ~doubly_linked_list_node() {}

        /**
         * @brief creates a node with value v
         */
        doubly_linked_list_node(T v) {
            this->v = v;
            pr = sc = NULL;
        };

        /**
         * @brief creates a node with value v, predecessor pr and successor sc
         * @param pr a node
         * @param sc a node
         */
        doubly_linked_list_node(T v, doubly_linked_list_node* pr, doubly_linked_list_node* sc) {
            this->v = v;
            this->pr = pr;
            this->sc = sc;
        }
    };

    protected:
    doubly_linked_list_node *hd = NULL; // first node
    doubly_linked_list_node *tl = NULL; // last node
    uint64_t s = 0; // size

    /**
     * @brief creates a copy of another list
     * @param other the list to copy
     */
    void copy_from_other(const doubly_linked_list& other) {
        if (other.s > 0) {
            this->s = other.s;
            doubly_linked_list_node* cur = other.hd;
            push_back(cur->v);

            while (cur->sc != NULL) {
                push_back(cur->v);
                cur = cur->sc;
            }
        }
    }

    /**
     * @brief moves a list into this list
     * @param other the list to move
     */
    void move_from_other(doubly_linked_list&& other) {
        this->hd = other.hd;
        this->tl = other.tl;
        this->s = other.s;

        other.hd = NULL;
        other.tl = NULL;
        other.s = NULL;
    }
    
    public:
    doubly_linked_list() = default;
    doubly_linked_list(doubly_linked_list&& other) {move_from_other(std::move(other));}
    doubly_linked_list(const doubly_linked_list& other) {copy_from_other(other);}
    doubly_linked_list& operator=(doubly_linked_list&& other) {move_from_other(std::move(other));return *this;}
    doubly_linked_list& operator=(const doubly_linked_list& other) {copy_from_other(other); return *this;}

    /**
     * @brief deletes all nodes of the list and the lsit
     */
    ~doubly_linked_list() {
        delete_nodes();
        hd = tl = NULL;
    }

    /**
     * @brief checks whether the list is empty
     * @return whether the list is empty
     */
    inline bool empty() {
        return s == 0;
    }

    /**
     * @brief returns the number of nodes in the list
     * @return number of nodes in the list
     */
    inline uint64_t size() {
        return s;
    }

    /**
     * @brief adjusts the size of the list
     * @param s new size
     */
    inline void set_size(uint64_t s) {
        this->s = s;
    }

    /**
     * @brief returns the head of the list
     * @return the head of the list
     */
    inline doubly_linked_list_node* head() {
        return hd;
    }

    /**
     * @brief adjusts the head of the list
     * @param n new head
     */
    inline void set_head(doubly_linked_list_node *n) {
        this->hd = n;
    }

    /**
     * @brief returns the tail of the list
     * @return the tail of the list
     */
    inline doubly_linked_list_node* tail() {
        return tl;
    }

    /**
     * @brief adjusts the head of the list
     * @param n new tail
     */
    inline void set_tail(doubly_linked_list_node *n) {
        this->tl = n;
    }

    /**
     * @brief inserts a node with value v before the head of the list
     * @param v value
     * @return the newly created node
     */
    inline doubly_linked_list_node* push_front(T &&v) {
        return push_front_node(new doubly_linked_list_node(v));
    }

    /**
     * @brief inserts a node with value v before the head of the list
     * @param v value
     * @return the newly created node
     */
    inline doubly_linked_list_node* push_front(T &v) {
        return push_front_node(new doubly_linked_list_node(std::move(v)));
    }

    /**
     * @brief inserts a node with value v after the head of the list
     * @param v value
     * @return the newly created node
     */
    inline doubly_linked_list_node* push_back(T &&v) {
        return push_back_node(new doubly_linked_list_node(v));
    }

    /**
     * @brief inserts a node with value v after the head of the list
     * @param v value
     * @return the newly created node
     */
    inline doubly_linked_list_node* push_back(T &v) {
        return push_back_node(new doubly_linked_list_node(std::move(v)));
    }

    /**
     * @brief inserts a node with value v before the node n2
     * @param v value
     * @param n2 a node in the list
     * @return the newly created node
     */
    inline doubly_linked_list_node* insert_before(T &&v, doubly_linked_list_node *n) {
        return insert_before_node(new doubly_linked_list_node(v),n);
    }

    /**
     * @brief inserts a node with value v before the node n2
     * @param v value
     * @param n2 a node in the list
     * @return the newly created node
     */
    inline doubly_linked_list_node* insert_before(T &v, doubly_linked_list_node *n) {
        return insert_before_node(new doubly_linked_list_node(std::move(v)),n);
    }

    /**
     * @brief inserts a node with value v after the node n2
     * @param v value
     * @param n2 a node in the list
     * @return the newly created node
     */
    inline doubly_linked_list_node* insert_after(T &&v, doubly_linked_list_node *n) {
        return insert_after_node(new doubly_linked_list_node(v),n);
    }

    /**
     * @brief inserts a node with value v after the node n2
     * @param v value
     * @param n2 a node in the list
     * @return the newly created node
     */
    inline doubly_linked_list_node* insert_after(T &v, doubly_linked_list_node *n) {
        return insert_after_node(new doubly_linked_list_node(std::move(v)),n);
    }

    /**
     * @brief inserts the node n1 before the node n2
     * @param n1 a node, that is not in the list
     * @param n2 a node in the list, n1 != n2
     */
    inline void insert_before_node(doubly_linked_list_node *n1, doubly_linked_list_node *n2) {
        if (n2 == hd) {
            hd = n1;
        } else {
            n2->pr->sc = n1;
            n1->pr = n2->pr;
        }
        n1->sc = n2;
        n2->pr = n1;
        s++;
    }

    /**
     * @brief inserts the node n1 after the node n2
     * @param n1 a node, that is not in the list
     * @param n2 a node in the list, n1 != n2
     */
    inline void insert_after_node(doubly_linked_list_node *n1, doubly_linked_list_node *n2) {
        if (n2 == tl) {
            tl = n1;
        } else {
            n2->sc->pr = n1;
            n1->sc = n2->sc;
        }
        n1->pr = n2;
        n2->sc = n1;
        s++;
    }

    /**
     * @brief inserts the node n before the head of the list
     * @param n a node, that is not in the list
     */
    inline void push_front_node(doubly_linked_list_node *n) {
        if (empty()) {
            hd = tl = n;
            s = 1;
        } else {
            insert_before_node(n,hd);
        }
    }

    /**
     * @brief inserts the node n after the tail of the list
     * @param n a node, that is not in the list
     */
    inline doubly_linked_list_node* push_back_node(doubly_linked_list_node *n) {
        if (empty()) {
            hd = tl = n;
            s = 1;
        } else {
            insert_after_node(n,tl);
        }
        return n;
    }

    /**
     * @brief concatenates the list l to the end of the list
     * @param l another list
     */
    inline void concat(doubly_linked_list<T> *l) {
        if (empty()) {
            if (!l->empty()) {
                hd = l->hd;
                tl = l->tl;
                s = l->s;
                l->disconnect_nodes();
            }
        } else if (!l->empty()) {
            tl->sc = l->hd;
            l->hd->pr = tl;
            tl = l->tail();
            s += l->s;
            l->disconnect_nodes();
        }
    }

    /**
     * @brief removes the node n from the list
     * @param n a node in the list
     */
    inline void remove_node(doubly_linked_list_node *n) {
        s--;
        if (n == hd) {
            hd = n->sc;
        } else if (n == tl) {
            tl = n->pr;
        }
        if (n->pr != NULL) {
            n->pr->sc = n->sc;
        }
        if (n->sc != NULL) {
            n->sc->pr = n->pr;
        }
    }

    /**
     * @brief disconnects the list from it's nodes
     */
    inline void disconnect_nodes() {
        hd = tl = NULL;
        s = 0;
    }

    /**
     * @brief deletes all nodes in the list
     */
    void delete_nodes() {
        if (!empty()) {
            doubly_linked_list_node *n = hd;
            for (uint64_t i=1; i<s; i++) {
                n = n->sc;
                delete n->pr;
            }
            hd = tl = NULL;
            s = 0;
        }
    }

    /**
     * @brief iterator for a list
     */
    class dll_it {
        protected:
        doubly_linked_list<T> *l; // the list, the iterator iterates through
        doubly_linked_list_node *cur; // the node the iterator points to

        public:
        /**
         * @brief creates an dl_it pointing to the node n in the list l
         * @param l a list
         * @param n a node in l
         */
        dll_it(doubly_linked_list<T> *l, doubly_linked_list_node *n) {
            this->l = l;
            this->cur = n;
        }

        /**
         * @brief deletes the iterator
         */
        ~dll_it() {
            l = NULL;
            cur = NULL;
        }

        /**
         * @brief checks whether the iterator can iterate forward
         * @return whether it can iterate forward
         */
        inline bool has_next() {
            return cur->sc != NULL;
        }

        /**
         * @brief checks whether the iterator can iterate backward
         * @return whether it can iterate backward
         */
        inline bool has_prev() {
            return cur->pr != NULL;
        }

        /**
         * @brief returns the value of the node, the iterator points to
         * @return the node, the iterator points to
         */
        inline doubly_linked_list_node* current() {
            return cur;
        }

        /**
         * @brief iterates forward, has_next() must return true
         * @return the node the iterator points to after iterating forward
         */
        inline doubly_linked_list_node* next() {
            cur = cur->sc;
            return cur;
        }

        /**
         * @brief iterates forward, has_pred() must return true
         * @return the node the iterator points to after iterating backward
         */
        inline doubly_linked_list_node* previous() {
            cur = cur->pr;
            return cur;
        }

        /**
         * @brief points the iterator to the node n
         * @param n a node in l
         */
        inline void set(doubly_linked_list_node *n) {
            cur = n;
        }
    };

    /**
     * @brief returns an iterator pointing to the node n
     * @param n a node in the list
     * @return an iterator
     */
    inline doubly_linked_list<T>::dll_it iterator(doubly_linked_list_node *n) {
        return doubly_linked_list<T>::dll_it(this,n);
    }

    /**
     * @brief returns an iterator pointing to the head of the list
     * @return an iterator
     */
    inline doubly_linked_list<T>::dll_it iterator() {
        return doubly_linked_list<T>::dll_it(this,hd);
    }
};