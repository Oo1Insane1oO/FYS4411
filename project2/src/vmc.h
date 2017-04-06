//////////////////////////////////////////////////////////////////////////////
// Header file for class vmc, see .cpp file for more information            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#ifndef VMC_H
#define VMC_H

#include "basis.h"

class VMC {
    private:
        Basis *b;
    
    public:
        VMC (Basis*);
        virtual ~VMC ();
};

#endif /* VMC_H */
