cdef extern from "<shared_mutex>" namespace "std" nogil:
    cdef cppclass shared_mutex:
        shared_mutex()
        void lock()
        void unlock()
        void lock_shared()
        void unlock_shared()