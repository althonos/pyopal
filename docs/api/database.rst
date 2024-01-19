Database
========

.. currentmodule:: pyopal


BaseDatabase
------------

.. autoclass:: pyopal.BaseDatabase
   :special-members: __init__, __len__, __getitem__
   :members:

   .. c:function:: size_t get_size(self)

      Return the number of elements in the database.

   .. c:function:: digit_t** get_sequences(self)

      Return a pointer to an array of sequence pointers.

   .. c:function:: int* get_lengths(self)
      
      Return a pointer to an array of lengths.


Database
--------

.. autoclass:: pyopal.Database(BaseDatabase)
   :special-members: __init__, __len__, __getitem__, __setitem__, __delitem__
   :members:
