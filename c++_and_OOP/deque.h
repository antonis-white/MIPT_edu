#include <iostream>
#include <vector>

//total snake_case
#define deque Deque

template <typename type>
class deque {
private:
    std::vector<type> buf;
    size_t size_ = 0;
    //position in buf where 0-sub-indexed element lies
    size_t shift = 0;
    //min/max sub-indexes in deque (may become negative)
    int beg_pos = 0;
    void realloc_up(size_t buf_size);
    //returns position in buffer by index
    size_t get_link(size_t index) const { return shift + index + beg_pos; }
    template <typename it_type>
    class iterator_ {
    private:
        //sub-index in deque (may become negative)
        int pos = 0;
        deque<type>* carrier = nullptr;
        friend class deque;
        iterator_(int pos, deque<type>* dq) : pos(pos), carrier(dq) {}
        //returns sub-index
    public:
        iterator_() = default;
        //returns index in deque
        size_t index() const { return pos - carrier->beg_pos; }
        iterator_& operator+=(int dlt) & { pos += dlt; return *this; }
        iterator_& operator-=(int dlt) & { pos -= dlt; return *this; }
        it_type& operator*() const { return carrier->buf[pos + carrier->shift]; }
        it_type* operator->() { return &(carrier->buf[pos + carrier->shift]); }
        void swap(iterator_& other) & { std::swap(pos, other.pos), std::swap(carrier, other.carrier); }
        template <typename it_type_v>
        bool operator==(const deque<type>::iterator_<it_type_v>& q) const {
            return pos == q.pos; 
        }
        template <typename it_type_v>
        bool operator<(const deque<type>::iterator_<it_type_v>& q) const {
            return pos < q.pos; 
        }
        template <typename it_type_v>
        bool operator>=(const deque<type>::iterator_<it_type_v>& q) const {
            return pos >= q.pos; 
        }
        template <typename it_type_v>
        bool operator<=(const deque<type>::iterator_<it_type_v>& q) const {
            return pos <= q.pos; 
        }
        template <typename it_type_v>
        bool operator>(const deque<type>::iterator_<it_type_v>& q) const {
            return pos > q.pos; 
        }
        template <typename it_type_v>
        bool operator!=(const deque<type>::iterator_<it_type_v>& q) const {
            return pos != q.pos; 
        }
        template <typename it_type_v>
        int operator-(const deque<type>::iterator_<it_type_v>& q) const {
            return pos - q.pos;
        }
        iterator_ operator+(int dlt) const {
            iterator_ res(*this);
            res += dlt;
            return res;
        }
        iterator_ operator-(int dlt) const {
            iterator_ res(*this);
            res -= dlt;
            return res;
        }
        friend iterator_ operator+(int dlt, const iterator_& it) {
            iterator_ res(it);
            res += dlt;
            return res;
        }
        friend iterator_ operator-(int dlt, const iterator_& it) {
            iterator_ res(it);
            res -= dlt;
            return res;
        }
        iterator_& operator++() & { *this += 1; return *this; }
        iterator_& operator--() & { *this -= 1; return *this; }
        iterator_ operator++(int) {
            iterator_ cp(*this);
            *this += 1;
            return cp;
        }
        iterator_ operator--(int) {
            iterator_ cp(*this);
            *this -= 1;
            return cp;
        }
        operator iterator_<const type>() { return iterator_<const type>(pos, carrier); }
        template <typename it_type_v>
        friend class iterator_;
    };
public:
    deque() = default;
    deque(size_t size_, const type& obj = type()) : buf(std::vector<type>(3 * size_, obj)), 
                                                    size_(size_), shift(size_), beg_pos(0) {}
    type& operator[](size_t index) { return buf[get_link(index)]; }
    const type& operator[](size_t index) const { return buf[get_link(index)]; }
    type& at(size_t index);            
    const type& at(size_t index) const;      
    size_t size() const { return size_; }
    void resize(size_t new_size, const type& obj = type());
    void reserve(size_t exp_size); 
    void push_back(const type& obj);       
    void push_front(const type& obj);     
    void pop_back() { --size_; }                     
    void pop_front() { --size_, ++beg_pos; }               
    using iterator = iterator_<type>;
    using const_iterator = iterator_<const type>;
    void insert(const const_iterator& it, const type& obj);
    void erase(const const_iterator& it); 
    iterator begin() { return iterator(beg_pos, this); }
    iterator end() { return iterator(beg_pos + size_, this); }
    const_iterator begin() const;
    const_iterator end() const;
    const_iterator cbegin() const;
    const_iterator cend() const;            
};

template <typename type>
typename deque<type>::const_iterator deque<type>::begin() const { 
    return const_iterator(beg_pos, const_cast<deque<type>*>(this)); 
}

template <typename type>
typename deque<type>::const_iterator deque<type>::end() const { 
    return const_iterator(beg_pos + size_, const_cast<deque<type>*>(this)); 
}

template <typename type>
typename deque<type>::const_iterator deque<type>::cbegin() const { 
    return const_iterator(beg_pos, const_cast<deque<type>*>(this)); 
}

template <typename type>
typename deque<type>::const_iterator deque<type>::cend() const { 
    return const_iterator(beg_pos + size_, const_cast<deque<type>*>(this)); 
}

template <typename type>
void deque<type>::realloc_up(size_t buf_size) {
    while ((buf_size % 3 != 0) || !buf_size) {
        ++buf_size;    
    }
    std::vector<type> nbuf(buf_size);
    size_t blank = buf_size / 3;
    auto nw_it = nbuf.begin() + blank;
    for (auto it = begin(); it != end(); ++it, ++nw_it) {
        *nw_it = *it;
    }
    shift = blank - beg_pos;
    buf.swap(nbuf);
}

template <typename type>
void deque<type>::insert(const const_iterator& it, const type& obj) {
    if (it == begin())
        return void(push_front(obj));
    if (get_link(size_) == buf.size()) 
        realloc_up(size_ * 3);
    int link = get_link(size_);
    for (size_t i = size_; i > it.index(); --i, --link) {
        std::swap(buf[link], buf[link - 1]);
    }
    ++size_;
    buf[link] = obj;
}

template <typename type>
void deque<type>::erase(const const_iterator& it) {
    if (it == begin()) 
        return void(pop_front());
    //where was mistake in testcase (erase(end()) is called in tests)
    if (it.index() >= size_) 
        return void(pop_back());        
    size_t link = get_link(it.index());
    --size_;
    for (size_t i = it.index(); i != size_; ++link, ++i) {
        std::swap(buf[link], buf[link + 1]);
    }
}

template <typename type>
type& deque<type>::at(size_t index) {
    if (index >= size_)
        throw std::out_of_range("index is out of deque");
    return buf[get_link(index)];
}

template <typename type>
const type& deque<type>::at(size_t index) const {
    if (index >= size_)
        throw std::out_of_range("index is out of deque");
    return buf[get_link(index)];
}

template <typename type>
void deque<type>::resize(size_t new_size, const type& obj) {
    if (new_size <= size_) 
        return void(size_ = new_size);
    if (get_link(beg_pos + new_size) > buf.size())
        realloc_up(3 * new_size);    
    size_t beg = get_link(size_), end = get_link(new_size);
    for (auto it = buf.begin() + beg; it != buf.begin() + end; ++it) {
        *it = obj;        
    }
    size_ = new_size;
}

template <typename type>
void deque<type>::reserve(size_t exp_size) {
    if (size_ < exp_size)
        realloc_up(exp_size);
}

template <typename type>
void deque<type>::push_back(const type& obj) {
    if (get_link(size_) == buf.size()) 
        realloc_up(3 * size_);
    buf[get_link(size_)] = obj;
    ++size_;
}

template <typename type>
void deque<type>::push_front(const type& obj) {
    if (get_link(0) == 0)
        realloc_up(3 * size_);
    --beg_pos;
    ++size_;
    buf[get_link(0)] = obj;
}

#undef deque