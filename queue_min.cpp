struct queueMin {

    queueMin() {}

    deque <int> q;

    int min() {
        return q.front();
    }

    void push(int x) {
        while (!q.empty() && q.back() > x)
            q.pop_back();
        q.push_back (x);
    }

    void pop(int x) {
        if (!q.empty() && q.front() == x)
            q.pop_front();
    }
};
