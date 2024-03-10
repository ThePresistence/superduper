#include <sys/category.h>
#include <sys/resource.h>
  /*   The following short routines provide a simple  
        interface for the Fortran user to some of the 
        system calls under Unicos. A synopsis of each 
        is given where appropriate. 
        A. Komornicki, Polyatomics Research Institute. 
                       Mountain View, CA.  August 1988. */ 
      MEMX(len)
      int *len;
   /*           Synopsis:  fl = memx(limit)  
      The value of limit passed may take on positive and
      negative values.  This will either increment or 
      decrement the current field length.   
      A value of zero will return the current field length. 
      Note:  Use of this routine is suggested only for cases 
      where a fixed heap is used at load time.   */ 
    { 
        return(sbreak(*len));
      }
      MEMLIM(lim) 
      long *lim; 
     /*   Synopsis:  fl = memlim(limit)  
                   or:  call memlim(limit) 
      This function or subroutine will return the current 
      process memory limit.  If invoked as a function, value 
      of input argument is ignored. 
      Note that if the system is not configured with limits 
      this routine will return a value of zero.  At this point 
      it is the user's responsibiity to use this information wisely. */ 
     {
        long  val;
        val = limit(C_PROC,0,L_MEM,-1);  
        if(val < 0)  {
         *lim = val;
         return(val);  
         }
         else {
           *lim = val*512;
           return(val*512);
          } 
        }
