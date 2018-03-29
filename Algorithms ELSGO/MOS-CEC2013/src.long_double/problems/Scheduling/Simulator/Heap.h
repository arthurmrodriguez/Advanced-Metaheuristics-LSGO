#ifndef _HEAP_H
#define _HEAP_H

#include <iostream>
#include <vector>
#include <algorithm>

//-----------------------------------------------------------------------------------------------------------------

template <typename T, class C> class Heap {

   public:

      typedef typename std::vector<T>::const_iterator const_iterator;

      Heap () {make_heap (mContainer.begin (), mContainer.end (), C ());}

      virtual ~Heap () {}

      const T& top  () const;
      void pop  ();
      void push (const T& x);

      void clear () {mContainer.clear ();}

      int  size  () const {return mContainer.size  ();}
      bool empty () const {return mContainer.empty ();}

      void sort () {sort_heap (mContainer.begin (), mContainer.end (), C ());}

      const_iterator begin () const {return mContainer.begin ();}
      const_iterator end   () const {return mContainer.end   ();}

      std::ostream& print (std::ostream& os) const;


   protected:

      std::vector<T> mContainer;

};


//-----------------------------------------------------------------------------------------------------------------

template <typename T, class C> const T& Heap<T,C>::top () const {

   if (!mContainer.size ()) {

      std::cerr << "Error: performing a top operation over an empty heap. Aborting..." << std::endl;
      abort ();

   }

   return mContainer.front ();

}

//-----------------------------------------------------------------------------------------------------------------

template <typename T, class C> void Heap<T,C>::pop () {

   if (!mContainer.size ()) {

      std::cerr << "Error: performing a pop operation over an empty heap. Aborting..." << std::endl;
      abort ();

   }

   pop_heap (mContainer.begin (), mContainer.end (), C ());
   mContainer.pop_back ();

}

//-----------------------------------------------------------------------------------------------------------------

template <typename T, class C> void Heap<T,C>::push (const T& x) {

   mContainer.push_back (x);

   push_heap (mContainer.begin (), mContainer.end (), C ());

}

//-----------------------------------------------------------------------------------------------------------------

template <typename T, class C> std::ostream& Heap<T,C>::print (std::ostream& os) const {

   const_iterator it;

   for (it = mContainer.begin (); it != mContainer.end (); it++)
      os << *it << " ";

   os << std::endl;

   return os;

}

//-----------------------------------------------------------------------------------------------------------------

#endif
