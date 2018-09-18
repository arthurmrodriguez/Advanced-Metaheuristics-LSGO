/* The following code declares classes to read from and write to
 * file descriptore or file handles.
 *
 * See
 *      http://www.josuttis.com/cppcode
 * for details and the latest version.
 *
 * - open:
 *      - integrating BUFSIZ on some systems?
 *      - optimized reading of multiple characters
 *      - stream for reading AND writing
 *      - i18n
 *
 * (C) Copyright Nicolai M. Josuttis 2001.
 * Permission to copy, use, modify, sell and distribute this software
 * is granted provided this copyright notice appears in all copies.
 * This software is provided "as is" without express or implied
 * warranty, and with no claim as to its suitability for any purpose.
 *
 * Version: Jul 28, 2002
 * History:
 *  Jul 28, 2002: bugfix memcpy() => memmove()
 *                fdinbuf::underflow(): cast for return statements
 *  Aug 05, 2001: first public version
 */
#ifndef BOOST_FDSTREAM_HPP
#define BOOST_FDSTREAM_HPP

#include <istream>
#include <ostream>
#include <streambuf>
// for EOF:
#include <cstdio>
// for memmove():
#include <cstring>

// low-level read and write functions
#ifdef _MSC_VER
# include <io.h>
#else
# include <unistd.h>
//extern "C" {
//    int write (int fd, const char* buf, int num);
//    int read (int fd, char* buf, int num);
//}
#endif


// BEGIN namespace BOOST
namespace boost {

/************************************************************
 * fdostream
 * - a stream that writes on a file descriptor
 ************************************************************/

class fdoutbuf : public std::streambuf {

   public:

      // Constructor
      fdoutbuf (int _fd) : fd (_fd) {}


   protected:

      int fd;    // File descriptor


      // Write one character
      virtual int_type overflow (int_type c) {

         if (c != EOF) {

            char z = c;

            if (write (fd, &z, 1) != 1)
               return EOF;

         }

         return c;

      }


      // Write multiple characters
      virtual std::streamsize xsputn (const char* s, std::streamsize num) {

         return write(fd,s,num);

      }

};


class fdostream : public std::ostream {

   public:

      fdostream (int fd) : std::ostream (0), buf (fd) {

         rdbuf (&buf);

      }


   protected:

      fdoutbuf buf;

};


/************************************************************
 * fdistream
 * - a stream that reads on a file descriptor
 ************************************************************/

class fdinbuf : public std::streambuf {

   public:

      /* constructor
       * - initialize file descriptor
       * - initialize empty data buffer
       * - no putback area
       * => force underflow()
       */

      fdinbuf (int _fd) : fd (_fd) {

         setg (buffer + pbSize,     // beginning of putback area
               buffer + pbSize,     // read position
               buffer + pbSize  );  // end position

      }


   protected:

      int fd;    // file descriptor

      /* data buffer:
       * - at most, pbSize characters in putback area plus
       * - at most, bufSize characters in ordinary read buffer
       */

      static const int pbSize  = 4;       // size of putback area
      static const int bufSize = 1024;    // size of the data buffer
      char buffer [bufSize + pbSize];     // data buffer

      // Insert new characters into the buffer
      virtual int_type underflow () {

#ifndef _MSC_VER
         using std::memmove;
#endif

         // Is read position before end of buffer?
         if (gptr () < egptr ())
            return traits_type::to_int_type (*gptr ());

         /* process size of putback area
          * - use number of characters read
          * - but at most size of putback area
          */

         int numPutback = gptr () - eback ();

         if (numPutback > pbSize)
            numPutback = pbSize;


         /* copy up to pbSize characters previously read into
          * the putback area
          */
         memmove (buffer + (pbSize - numPutback), gptr () - numPutback, numPutback);

         // Read at most bufSize new characters
         int num = read (fd, buffer + pbSize, bufSize);

         if (num <= 0) // ERROR or EOF
            return EOF;

         // Reset buffer pointers
         setg (buffer + pbSize - numPutback,   // beginning of putback area
               buffer + pbSize,                // read position
               buffer + pbSize + num);         // end of buffer

         // Return next character
         return traits_type::to_int_type (*gptr ());

      }

};


class fdistream : public std::istream {

   public:

      fdistream (int fd) : std::istream (0), buf (fd) {

         rdbuf (&buf);

      }


   protected:

      fdinbuf buf;

};

} // END namespace boost

#endif /*BOOST_FDSTREAM_HPP*/
