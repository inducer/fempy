/* $Id: triangmodule.c,v 1.12 2003/07/21 14:37:27 pletzer Exp $ */

/*  This code interfaces ellipt2d directly with "triangle", a general        
    purpose triangulation code. In doing so, Python data structures          
    are passed to C for  use by "triangle" and the results are then          
    converted back. This was accomplished using the Python/C API.            
                                                                           
    I. Lists of points,edges,holes, and regions must normally be 
     created for proper triangulation. However, if there are no holes,
     the edge list and hole list need not be specified (however, a 
     measure of control over the shape of the convex hull is lost; 
     triangle decides). In additon an intial area constraint for          
     triangles must be provided upon initialization of an Ireg2tri           
     object.

     points list: [(x1,y1),(x2,y2),...] (Tuples of doubles)                  
     edge list: [(Ea,Eb),(Ef,Eg),...] (Tuples of integers)                   
     hole list: [(x1,y1),...](Tuples of doubles, one inside each hole region)
     regionlist: [ (x1,y1,index,area),...] (Tuple of doubles)                

     This is all that is needed to accomplish an intial triangulation.       
     A dictionary in the form of the Node class data structure is            
     returned and passed to a Node object.                                   
                                                        
          A. To deal with Dirichlet boundary conditions a method           
     setUpDirchlet(...) is supplied. This takes a Dirichlet boundary       
     object and sets the second of four point attributes in triangle       
     to the boundary value. The first attribute is reserved for accounting 
     of unwanted interpolation of boudary values by triangle(for           
     each boundary node requiring Dirichlet boundary conditions it is      
     set to 1). After triangulation and/or refinement the method           
     updateDirichlet(...) is used to update the original dirichlet         
     Boundary object in the driver. It only returns points with            
     first attribute equal to 1 and on the boundary. Thus, using triangle  
     to deal with boudary value alterations due to added points during     
     refinement requires that BC's be specified only on a boundary and     
     not the interior of a domain. Otherwise,an update will ommit          
     non-boundary points. These points should be set afterwards.           
     (Triangle forms a smooth interpolation of sharp changes in boundary     
     values when intermediate points are added)                           
          B. Robbins (and its subset Neumann) boundary condtions are not      
     treated like the Dirichlet type(Mainly because Triangle lacks an      
     attribute list for edges). Users must employ the setAttributes(...)   
     and getAttributes(...) methods. These are specified using the         
     third and fourth point attributes. They may be used as seen fit.      
                                                                           
    II.Refinement of a mesh is accomplished by storing the original          
    output in a global variable(struct triangulateio out;). When the         
    refinement routine is called, a new area parameter is passed to the      
    triagulaton routine, as well as the 'r' switch to denote refinement.     
    A new mesh is returned in the triagulateio struct next, then the         
    memory for the old output is released and out is made to point to        
    the data of next.                                                        
                    
    Since out is global, the user can be dealing with only one mesh at a     
    time. However, multiple meshes can be used if  everything that needs     
    to be accomplished with a previous mesh is already done. Out will        
    simply point to the new mesh.                                            


    Compilation can be done as follows                                       
                    
    c++ -I/usr/local/include/python1.5 -shared -o triangmodule.so            
    triangmodule.c triangle.o                                                */


#ifdef SINGLE  
#define REAL float 
#else  
#define REAL double 
#endif  

#include "triangle.h" 


  static PyMethodDef triang_methods[] = {
    {(char *)"triangulate", triang_triangulate, 1 ,(char *)"Uses triangle to generate an unstructured mesh.(<PyString area>,<PyString mode>)-> Node object's data strucutre."},
    {(char *)"refine", triang_refine, 1, (char *)"Refines a previously generated mesh.(<PyString area>,<PyString mode>)-> Node object's data structure."},
    {(char *)"setUpConvexHull", triang_setUpConvexHull, 1, (char *)"Passes hull points to triangle.All arguments are lists.(<hull>,<seglist>,<holelist>,<regionlist>)-><hull>"},
    {(char *)"setUpDirichlet", triang_setUpDirichlet, 1, (char *)"Uses dirchlet boundary object to set point attributes in triangle.(<dB dictionary data strucuture>)-><dB>"},
    {(char *)"updateDirichlet", triang_updateDirichlet, 1, (char *)"Returns updated boundary value information(post triangulation),()-><dB dictionary data structure>"},
    {(char *)"setAttributes", triang_setAttributes, 1, (char *)"Directly sets point attributes,(<Dict of point attribtes>)-><Dictionary of point attributes>" },
    {(char *)"getAttributes", triang_getAttributes, 1, (char *)"Gets point attributes,()-><Dictionary of point attributes>"},
    {(char *)"delete", triang_delete, 1, (char *)"Frees memory allocated by triangle,()->Py_None"},
    {NULL,NULL}
  };

  void inittriang(){
    Py_InitModule((char *)"triang",triang_methods);
  }














































































