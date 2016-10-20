#include <vmdstream/vmdstream.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

//draw a bounding box so that the frame fits in the same area
//otherwise it will resize depending on where the particles are
/*void draw_bounding_box(
double xmin, double xmax, 
double ymin, double ymax, 
double zmin, double zmax)
{
    ofstream os("bounding_box.xyz");
    os << "10\n\n";
    os << "C\t" << xmin << "\t" << ymin << "\t" << zmin << "\n";
    os << "C\t" << xmin << "\t" << ymax << "\t" << zmin << "\n";
    os << "C\t" << xmax << "\t" << ymax << "\t" << zmin << "\n";
    os << "C\t" << xmax << "\t" << ymin << "\t" << zmin << "\n";
    os << "C\t" << xmin << "\t" << ymin << "\t" << zmin << "\n";
    os << "C\t" << xmin << "\t" << ymin << "\t" << zmax << "\n";
    os << "C\t" << xmin << "\t" << ymax << "\t" << zmax << "\n";
    os << "C\t" << xmax << "\t" << ymax << "\t" << zmax << "\n";
    os << "C\t" << xmax << "\t" << ymin << "\t" << zmax << "\n";
    os << "C\t" << xmin << "\t" << ymin << "\t" << zmax << "\n";
    os.close();     
}*/

#define MYPI 3.14159265

// box info
struct TriclinicBox{
    double a;
    double b;
    double c;
    double alpha;
    double beta;
    double gamma;

    TriclinicBox(double a1x, double a2x, double a2y, double a3x, double a3y, double a3z){
        a = sqrt(a1x*a1x);
        b = sqrt(a2x*a2x +a2y*a2y);
        c = sqrt(a3x*a3x + a3y*a3y + a3z*a3z);

        alpha = (180./MYPI)*acos((a2x*a3x + a2y*a3y)/(b*c));
        beta =  (180./MYPI)*acos(a3x/c);
        gamma = (180./MYPI)*acos(a2x/b);

        std::cout << "Volume of box is: " << getVolume() << endl;
        std::cout << a << "\t" << b << "\t" << "\t" << c << "\n"
                  << alpha << "\t" << beta << "\t" << gamma << "\n";

    }

    TriclinicBox(const char* filename="box_config.dat"){
        double a1x, a2x,a2y, a3x, a3y, a3z;
        ifstream myfile(filename);
        if (myfile.good()){
            myfile >> a1x >> a2x >> a2y >> a3x >> a3y >> a3z;
            //implememnt try and catch block

        }
        else{
            std::cerr << "Unable to read box config\n";
            exit(-2);
        }
        *this = TriclinicBox(a1x,a2x,a2y,a3x,a3y,a3z);
    }

    double getVolume(){
        double ca = cos(alpha*MYPI/180.);
        double cb = cos(beta*MYPI/180.);
        double cg = cos(gamma*MYPI/180.);
        return a*b*c*sqrt(1-ca*ca - cb*cb- cg*cg + 2*ca*cb*cg);
    }
};

int main (int argc, char** argv)
{        
    //initialize vmd; assumes that typing "vmd" in the terminal launches a VMD window
    int port = 5500;
      vmdsock_t vmdsock = newvmdsock("vmd", port);
    vmdstream vmdscript(vmdsock);

    //set up the canvas
    /*vmdscript << "axes location off\n";
    vmdscript << "color scale method BGR\n";
    vmdscript << "menu main off\n";
    vmdscript << "display resize 800 800\n"; 

    int nstep = 50;
    int pad = 10000;

    //initialize the position of the particle
    double x = 0.0, y = 0.0, z = 0.0;

    //draw the bounding box of size 4x4x4 centered at the origin
    draw_bounding_box(-2, 2, -2, 2, -2, 2);*/
    int nstep = 1;
    int pad = 10000;
    std::string myfile("lattice.xyz");
    std::ifstream vmf(myfile.c_str());
    int nmol=0;
    if (vmf.good()){
        std::string str;
        getline(vmf,str);
        nmol = std::stoi(str);
        nmol /= 4;
        vmf.close();

    }
    else exit(-1);

    TriclinicBox tBox("config.dat");

    for (int t=0; t<nstep; t++) {

        vmdscript << "draw delete all" << endl;          
        vmdscript << "color Display Background iceblue" << endl;
        vmdscript << "draw material Opaque" << endl;

        //normally one would load in data from a file or update a simulation here
        //I will simply do a random walk for a single sphere
        /*x += (2.0*rand()/(double)RAND_MAX-1.0)/10.0;
        y += (2.0*rand()/(double)RAND_MAX-1.0)/10.0;
        z += (2.0*rand()/(double)RAND_MAX-1.0)/10.0;*/

        // load file
        vmdscript << "mol addrep 0" << endl;
        vmdscript << "mol new {"<<myfile<<"} type {xyz} first 0 last -1 step 1 waitfor 1" << endl;


        // set up config
        vmdscript << "topo clearbonds" << endl;

        //set the properties of the shape
        int color = rand()%33;
        double radius = 0.25;


        double ax, bx, cx, alpha, beta, gamma;
        ax = tBox.a;
        bx = tBox.b;
        cx = tBox.c;
        alpha=tBox.alpha;
        beta=tBox.beta;
        gamma=tBox.gamma;
        string spacer=", ";

        std::cout << ax << "\t" << to_string(ax) << "\t" << tBox.a << endl;

        //draw the shape (do this in a for loop for all the shapes in the system)
        //vmdscript << "draw color " << color << endl;
        //vmdscript << "menu graphics on" << endl;
        vmdscript << "mol modstyle 0 1 CPK 1.000000 0.300000 12.000000 12.000000" << endl;
        vmdscript << "mol modstyle 0 1 CPK 0.600000 0.500000 15.000000 15.000000" << endl;
        //vmdscript << "pbc set { 26.270269999999996, 26.933820476852251, 24.308296665720896, "
        //    << "89.209926761906445, 83.418417942993528, 83.418417942993528 } -first now -all " << endl;
        string set_box = "pbc set { " + to_string(ax) + spacer + to_string(bx) + spacer + to_string(cx) +
                          spacer + to_string(alpha) + spacer + to_string(beta) + spacer + to_string(gamma) +
                          " } -first now -all ";
        vmdscript << set_box << endl;
        std::cout << set_box << endl;
        //vmdscript << "menu graphics off" << endl;

        for (unsigned ind =0; ind < 3; ind++ ){
            vmdscript << "for { set x 0 } { $x < "<< nmol <<" } { incr x }"
                  << " { puts [ topo addbond [ expr 4*$x + " << ind << "]"
                  << " [ expr 4*$x + " << ind+1 << " ] ] }" << endl;
            if (ind < 2){
                vmdscript << "for { set x 0 } { $x < "<< nmol <<" } { incr x }"
                  << " { puts [ topo addangle [ expr 4*$x + " << ind << "]"
                  << " [ expr 4*$x + " << ind+1 << " ] [ expr 4*$x + "
                  << ind+2 << " ] ] }" << endl;
                  if (ind < 1){
                    vmdscript << "for { set x 0 } { $x < "<< nmol <<" } { incr x }"
                    << " { puts [ topo adddihedral [ expr 4*$x + " << ind << "]"
                    << " [ expr 4*$x + " << ind+1 << " ] [ expr 4*$x + "
                    << ind+2  << " ] [ expr 4*$x + " << ind+3 << " ] ] }" << endl;
                  }
            }
        }
      


        //set the properties of the shape
        /*int color = rand()%33;
        double radius = 0.25;

        //draw the shape (do this in a for loop for all the shapes in the system)
        vmdscript << "draw color " << color << endl;
        vmdscript << "draw sphere {" << x << " "<< y << " "<< z << "} radius " << radius << endl;

        //now draw the bounding box
        vmdscript << "mol new bounding_box.xyz" << endl;
        vmdscript << "mol modcolor 0 top ColorID 8" << endl;
        vmdscript << "mol modstyle 0 top DynamicBonds 4.10000 0.0500000 6.000000" << endl;
        vmdscript << "scale by 0.8" << endl;
        vmdscript << "rotate y by " << t*1 << endl;
        vmdscript << "rotate x by 5" << endl;*/

//#ifdef DNDEBUG
  //      std::cin.get() ;
//#endif

        //render the frame; pad variable ensures the image files "ls" in the proper order
        vmdscript << "render snapshot frame" << pad + t << ".tga" << endl;

        //write lammpsdata file
        vmdscript << "topo writelammpsdata large_"<<t<<".lmp molecular" << endl;


        //delete everything and start a new frame
        vmdscript << "mol delete all" << endl;
    }

    //quit vmd and exit
    vmdscript << "quit" << endl;        
    vmdscript.flush();
    closevmdsock(vmdsock);

    //try to post process file
    system("./postprocess.sh");


    exit(1);
}