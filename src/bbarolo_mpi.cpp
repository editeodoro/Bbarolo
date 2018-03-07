
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <bbarolo.hh>
#include <Arrays/param.hh>

int main (int argc, char *argv[]) {

    if (argc!=3) {
        std::cerr << "Usage: BBarolo_MPI -l listfile \n";
        return EXIT_FAILURE;
    }

    if (std::string(argv[1])!="-l") {
        std::cerr << "Usage: BBarolo_MPI -l listfile \n";
        return EXIT_FAILURE;
    }
    
    int rank=0, nprocs=1;
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string listfile = std::string(argv[2]);
    std::vector<std::string> parfiles;

    // Reading parameter list file
    std::ifstream file(listfile.c_str());
    if(!file) std::cerr << "File " << listfile << " does not exist! \n";
    std::string s;
    while(file.good()) {
        getline(file, s);
        if (s!="") parfiles.push_back(s);
    }
    file.close();
               
    // Dividing parameter files among MPI processes
    int lsize = parfiles.size()/nprocs;
    int res_lsize = parfiles.size()%nprocs;
        
    int start[nprocs+1];
    for (int i=0; i<nprocs; i++) {
        if (i<res_lsize) start[i]=(lsize+1)*i;
        else start[i]=i*lsize+res_lsize;
    }
    start[nprocs]=parfiles.size();
    
    // Redirecting std::cout to files
    std::filebuf buf;
    std::string pout = "pout."+to_string(rank);
    buf.open(pout, std::ios::out );
    auto oldbuf = std::cout.rdbuf(&buf);
    
    // Main loop over parameter files
    for (auto i=start[rank]; i<start[rank+1]; i++) {
       
        Param *par = new Param;
        par->readParams(parfiles[i]);
        std::cout << *par;
        
        if (!BBcore(par)) {
            if(par->getListSize()-i>1) std::cout << "Skipping to next file...\n";
            else {std::cout << "Exiting ...\n\n"; return EXIT_FAILURE;}
            delete par;
            continue;
        }

        delete par;
    }

    // Back to std::cout
    std::cout.rdbuf(oldbuf);
	MPI_Finalize();
    
}
