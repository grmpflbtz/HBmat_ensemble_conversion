#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

int main(int argc, char *argv[]) {

    if( argc !=6 ) {
        std::cout << "input format:" << std::endl << "ln g(U) input_file" << std::endl << "HBmatrix input_file" << std::endl << "minimum_temperature" << std::endl << "maximum_temperature" << std::endl << "temperature_step" << std::endl;
        return -10;
    }
    std::cout << std::setprecision(15) 
    << "parameters:" << std::endl 
    << "lng_file:      " << argv[1] << std::endl 
    << "HBmatrix_file: " << argv[2] << std::endl 
    << "Tmin:    " << argv[3] << std::endl 
    << "Tmax:    " << argv[4] << std::endl 
    << "Tstep:   " << argv[5] << std::endl;


    double Tmin=-500, Tmax=500, Tstep=1;
    int NumTSteps;
    int i;      // itteration variable

    Tmin = atof(argv[3]);
    Tmax = atof(argv[4]);
    Tstep = atof(argv[5]);
    NumTSteps = (Tmax-Tmin)/(Tstep);

    std::cout << "#Tsteps: " << NumTSteps << std::endl;

    // reading lng(U) file
    std::ifstream inpDOS1(argv[1]);
    std::string first_line;
    double x,y,z;
    double EnMin=1e30, EnMax=-1e30;
    int Ecount=0;

    if(!inpDOS1.is_open()) {
        std::cout << " --- ERROR --- missing ln g(U) input_file:"  << argv[1] << std::endl;
        return -10;
    }
    getline(inpDOS1, first_line);
    for( i=0; !(inpDOS1.eof()); i++) {
        inpDOS1 >> x >> x >> y >> z >> y;
        if(x<EnMin) {EnMin=x;}
        if(x>EnMax) {EnMax=x;}
        Ecount++;
    }
    Ecount--;
    std::cout << "Emin:    " << EnMin << std::endl 
    << "Emax:    " <<EnMax << std::endl 
    << "Ecount:  " << Ecount << std::endl;
    inpDOS1.close();

    double *En, *lng;

    En = new double [Ecount];
    lng = new double [Ecount];

    std::ifstream inpDos2(argv[1]);
    getline(inpDos2, first_line);
    inpDos2 >> x >> x >> y >> z >> y;
    for( i=0; !(inpDos2.eof()); i++ ) {
        En[i] = x;
        lng[i] = z;
        inpDos2 >> x >> x >> y >> z >> y;
    }
    inpDos2.close();

    /*for( i=0; i<Ecount; i++) {
        std::cout << En[i] << " " << lng[i] << std::endl;
    }
    */

    // reading HBmatrix file
    std::ifstream inpHBM(argv[2]);
    if(!inpHBM.is_open()) {
        std::cout << " --- ERROR --- missing HBmatrix input_file:"  << argv[2] << std::endl;
        delete [] En;
        delete [] lng;
        return -10;
    }
    int AAnum, nrow;
    char delim = ' ';
    std::string line, line_seg;
//    std::stringstream ss_line;

    while(true) {
        getline(inpHBM, first_line);
        if(first_line[0] == '#') {
            continue;
        }
        else{
            AAnum = std::stoi(first_line);
            break;
        }
    }
    
    
    //inpHBM >> AAnum;

    if( AAnum == 0 ) {
        std::cout << " --- ERROR --- no. of AA = 0" << std::endl
                  << "     ---->  probably AAnum missing before matrices" << std::endl;
        delete [] En;
        delete [] lng;
        return -10;
    }

    double **HBmatE;
    double **HBmatT;
    double *HBline;
    HBmatE = new double*[AAnum*Ecount];
    for( i=0; i<AAnum*Ecount; i++ ) { HBmatE[i] = new double [AAnum]; }
    HBmatT = new double*[AAnum*NumTSteps];
    for( i=0; i<AAnum*NumTSteps; i++ ) { HBmatT[i] = new double [AAnum]; }
    HBline = new double[AAnum];

    nrow=0;
    std::getline(inpHBM, line);
    while( std::getline(inpHBM, line) ) {

        //std::cout << line <<  "    ->    ";

        std::stringstream ss_line(line) ;
        for( i=0; i<AAnum; i++) {
            getline(ss_line, line_seg, delim);


            //std::cout << line_seg << " ";


            if( line_seg.compare("-nan") == 0 || line_seg.compare("nan") == 0 ) {
                x = -1;
            }
            else {
                x = std::stod(line_seg, NULL);
            }
            HBmatE[nrow][i] = x;
        }



        //std::cout << std::endl;



        nrow++;
    }
    inpHBM.close();

/*
    for( int i=0; i<AAnum*AAnum*Ecount; i++ ) {
        std::cout << HBmatE[i/AAnum][i%AAnum] << " ";
        if((i+1)%AAnum == 0) {
            std::cout << std::endl;
        }
    }
*/


    double *Zpart, *HBTpart;
    double Z, ZpartMax, ZpartMax1;
    double Tcur;
    int maxj = Ecount;

    Zpart = new double [Ecount];
    HBTpart = new double [Ecount];

    for( i=0,Tcur=Tmin; Tcur<Tmax; i++,Tcur+=Tstep ) {

        //std::cout << "Tcur = " << Tcur << std::endl;

        // calculate Z(Tcur)
        ZpartMax = -1e30;
        for( int j=0; j<maxj; j++ ) {
            Zpart[j] = lng[j] - En[j]/Tcur;
            if( Zpart[j] > ZpartMax ) { ZpartMax = Zpart[j]; }
        }
        for( int j=0; j<maxj; j++ ) {
            Zpart[j] -= ZpartMax;

            //std::cout << j << "\t" << Zpart[j] << std::endl;

        }
        ZpartMax1 = ZpartMax;

        std::sort(Zpart, Zpart+maxj);
        Z=0;
        for( int j=0; j<maxj; j++ ) {
            Z += exp( Zpart[j] );

            //std::cout << j << "\t" << Z << std::endl;

        }

        // calculate HBTpart -> HBmaxT

        for( int m=0; m<AAnum; m++ ) {
            for( int n=0; n<AAnum; n++ ) {
                for( int j=0; j<maxj; j++) {
                    if( HBmatE[j*AAnum+m][n] == -1 ) { 
                        HBTpart[j] = 0; 
                    } else { 
                        HBTpart[j] = HBmatE[j*AAnum+m][n] * exp(lng[j]-En[j]/Tcur - ZpartMax1);
                    }
                }
                std::sort(HBTpart, HBTpart+maxj);
                HBmatT[i*AAnum+m][n] = 0;
                for( int j=0; j<maxj; j++ ) {

                    //if(i>0){ 
                    //    std::cout << HBTpart[j] << std::flush; 
                    //}

                    HBmatT[i*AAnum+m][n] += HBTpart[j];
                }
                HBmatT[i*AAnum+m][n] /= Z;
            }
        }
    }

    std::ofstream outpHBM("HBmat_P20_T.dat");
    outpHBM << "# hydrogen bond matrices in the canonical ensemble" << std::endl 
            << "# Tmin=" << Tmin << ", Tmax=" << Tmax << ", dT=" << Tstep << std::endl
            << AAnum << std::endl;
    for( i=0; i<NumTSteps; i++ ) {
        for( int m=0; m<AAnum; m++ ) {
            for( int n=0; n<AAnum; n++ ) {
                outpHBM << HBmatT[i*AAnum+m][n] << " ";
            }
            outpHBM << std::endl;
        }
    }
    outpHBM.close();


    delete [] En;
    delete [] lng;
    for( i=0; i<AAnum*Ecount; i++ ) { delete [] HBmatE[i]; }
    delete [] HBmatE;
    for( i=0; i<AAnum*NumTSteps; i++ ) { delete [] HBmatT[i]; }
    delete [] HBmatT;
    delete [] HBline;
    delete [] Zpart;
    delete [] HBTpart;

    return 0;
}