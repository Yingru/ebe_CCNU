#include"CEos.h"

CEOSL::CEOSL(){
}

void CEOSL::init()
{
    for(int i=0; i<4; i++){
        char EPname[64], Tname[64];

        sprintf(EPname,"s95p-v1/s95p-v1_dens%d.dat",i+1);
        sprintf( Tname,"s95p-v1/s95p-v1_par%d.dat", i+1);
        ifstream finep(EPname);
        ifstream fint(Tname);
        double e,p,s,T,buff;

        finep>>t[i].e0>>t[i].de>>t[i].Ne;
        fint>>buff>>buff>>buff;

        for(int n=0; n<t[i].Ne; ++n){
            finep>>e>>p>>s>>buff>>buff;
            fint>>T>>buff>>buff;
            t[i].e.push_front(e);
            t[i].p.push_front(p);
            t[i].s.push_front(s);
            t[i].T.push_front(T);
        }
    }//Construct the 4 tables for e, p, T
}

double CEOSL::lin_int(CD &x1, CD &x2, CD &y1, CD &y2, CD &x){
    double r= (x-x1)/(x2-x1);
    return y1*(1-r) + y2*r;
}

//////////////////////////////////////////////////////////////////////
double CEOSL::P(CD &eps){
    deque<double>::size_type i;
    double Pr;
    double emin = t[0].e0;
    double emax = t[3].e0+t[3].de*(t[3].Ne-1.0);
    if(eps>= emin && eps<= emax){
        for(int j=3; j>=0; j--){
            if( eps>=t[j].e0 ){
                i = unsigned( (eps-t[j].e0)/t[j].de );
                Pr = lin_int(t[j].e[i], t[j].e[i+1], t[j].p[i], t[j].p[i+1], eps);
                break;
            }
        }// Look for table position and do interpolation
    }
    else if(eps<emin)Pr= lin_int(0.0, t[0].e0, 0.0, t[0].p[0], eps);

    else if(eps>emax){
        double e1=t[3].e0+t[3].de*(t[3].Ne-2.0);
        double e2=t[3].e0+t[3].de*(t[3].Ne-1.0);
        double p1=t[3].p[t[3].Ne-2];
        double p2=t[3].p[t[3].Ne-1];
        Pr= lin_int(e1, e2, p1, p2, eps);
        //Pr=0.3327*eps-0.3223*pow(eps,0.4585)-0.003906*eps*exp(-0.05697*eps);
    }
    return Pr;
}
//////////////////////////////////////////////////////////////////////
double CEOSL::T(CD &eps){
    deque<double>::size_type i;
    double T;
    double emin = t[0].e0;
    double emax = t[3].e0+t[3].de*(t[3].Ne-1.0);
    if(eps>= emin && eps<= emax){
        for(int j=3; j>=0; j--){
            if( eps>=t[j].e0 ){
                i = unsigned( (eps-t[j].e0)/t[j].de );
                T = lin_int(t[j].e[i], t[j].e[i+1], t[j].T[i], t[j].T[i+1], eps);
                break;
            }
        }// Look for table position and do interpolation
    }
    else if(eps<emin)T = lin_int(0.0, t[0].e0, 0.0, t[0].T[0], eps);
    else if(eps>emax){
        double e1=t[3].e0+t[3].de*(t[3].Ne-2.0);
        double e2=t[3].e0+t[3].de*(t[3].Ne-1.0);
        double T1=t[3].T[t[3].Ne-2];
        double T2=t[3].T[t[3].Ne-1];
        T= lin_int(e1, e2, T1, T2, eps);
        //Pr=0.3327*eps-0.3223*pow(eps,0.4585)-0.003906*eps*exp(-0.05697*eps);
    }
    return T;
}
//////////////////////////////////////////////////////////////////////
double CEOSL::S(CD &eps){
    deque<double>::size_type i;
    double s;

    double emin = t[0].e0;
    double emax = t[3].e0+t[3].de*(t[3].Ne-1.0);
    if(eps>= emin && eps<= emax){
        for(int j=3; j>=0; j--){
            if( eps>=t[j].e0 ){
                i = unsigned( (eps-t[j].e0)/t[j].de );
                s = lin_int(t[j].e[i], t[j].e[i+1], t[j].s[i], t[j].s[i+1], eps);
                break;
            }
        }// Look for table position and do interpolation
    }
    else if(eps<emin) s = lin_int(0.0, t[0].e0, 0.0, t[0].s[0], eps);

    else if(eps>emax){
        double e1=t[3].e0+t[3].de*(t[3].Ne-2.0);
        double e2=t[3].e0+t[3].de*(t[3].Ne-1.0);
        double s1=t[3].s[t[3].Ne-2];
        double s2=t[3].s[t[3].Ne-1];
        s = lin_int(e1, e2, s1, s2, eps);
        //Pr=0.3327*eps-0.3223*pow(eps,0.4585)-0.003906*eps*exp(-0.05697*eps);
    }
    return s;

}
//////////////////////////////////////////////////////////////////////

/* EOSL https://wiki.bnl.gov/TECHQM/index.php/QCD_Equation_of_State */

CEos::CEos()
{
    IEOS=0;  //default construct function
}

CEos::CEos(int ieos):IEOS(ieos)
{
    //construct function
    if(ieos==EOSL2)eosl.init();
}

void CEos::SetEos(int ieos)
{
    IEOS=ieos;
    if(ieos==EOSL2)eosl.init();
}

double CEos::P(CD &eps, CD &bd) //unit Gev/fm^3
{   
    double Pr;
    double ee;
    double Aperp;
    double Cperp;

    CD eosl_e1=0.5028563305441270;
    CD eosl_e2=1.62;
    CD eosl_e3=1.86;
    CD eosl_e4=9.9878355786273545;

    switch(IEOS){
        case EOSI:        
            Pr=eps/3.0;
            break;
        case EOSQ:	
            break;
        case SMEOSQ://close to peter sm-eosq   
            ee=eps/hbarc;     //unit fm^-4
            if(ee>30.0)Pr=hbarc*(0.334*ee-2.46085238);
            else if(ee>5.0)
                Pr=hbarc*(-2.46085238+0.167*log(2.696674E7+exp(2.0*ee)));
            else
                Pr=hbarc*max(0.0, 0.152000165*ee				\
                        -0.0152*log(4.19314E10+exp(10.0*ee))		\
                        -0.0248909*exp(-2.7*ee) +0.3966722716141761);
            break;
        case EOSL://Lattice EOS               
            if(eps<eosl_e1)Pr=0.3299*(exp(0.4346*eps)-1.0);
            else if(eps>=eosl_e1 && eps<eosl_e2)Pr=1.024*1.0E-7*exp(6.041*eps)+0.007273+0.14578*eps;
            else if(eps>=eosl_e2 && eps<eosl_e3)Pr=0.30195*exp(0.31308*eps)-0.256232;
            else if(eps>=eosl_e3 && eps<eosl_e4)Pr=0.332*eps-0.3223*pow(eps,0.4585)-0.003906*eps*exp(-0.05697*eps)\
                +0.1167*pow(eps,-1.233)+0.1436*eps*exp(-0.913*eps);
            else if(eps>eosl_e4)Pr=0.3327*eps-0.3223*pow(eps,0.4585)-0.003906*eps*exp(-0.05697*eps);
            //	    ee=eps/hbarc;     //unit fm^-4
            //	    Aperp=1.0/(exp((4.0-ee)/1.0)+1.0);   
            //	    Cperp=1.0/(exp((1.8-ee)/1.5)+1.0);
            //	    Pr=hbarc*max(0.0,(-2.2+2.655*log(1.18136+exp(0.1111*ee)))*Aperp	 \
            //		+0.18*ee*(1.0-Cperp)+(0.40*Cperp)*(1-Aperp)	 \
            //		-0.0886002645519078);
            break;

        case EOSL2://Lattice EOS s95p-v1
            Pr=eosl.P(eps);
            break;

        default:
            Pr=eps/3.0;
    }
    return Pr;
}

double CEos::T(CD &eps, CD &bd)
{	double Temp;
    double gt=169.0/4.0;         //Total freedom of Quarks and Gluon
    double a=(gt*pi*pi)/90.0;
    double ConsT=pow(1.0/(3.0*a),0.25);
    double ee,Aper, Bper, Cper;
    switch(IEOS){
        case EOSI:
            ee=eps/hbarc;
            Temp=hbarc*ConsT*pow(ee,0.25);
            break;
        case EOSQ:
            break;
        case SMEOSQ://Peter EOS SM-EOSQ
            ee=eps/hbarc;
            Bper=1.0/(exp((8.0-ee)/3.0)+1.0);
            Cper=1.0/(exp((5.0-ee)/15.0)+1.0);
            Temp=hbarc*(0.0+pow(0.07*ee,0.25)*Bper	\
                    +pow(6.2*ee,0.16)*(1-Cper)*(1-Bper));
            break;
        case EOSL:
            if(eps<0.5143939846236409)Temp=0.203054*pow(eps,0.30679);
            else Temp=(eps+P(eps,bd))/S(eps,bd);
            //    ee=eps/hbarc;
            //    Aper=1.0/(exp((60.0-ee)/60.0)+1.0);   
            //    Bper=1.0/(exp((1.5-ee)/7.0)+1.0);
            //    Cper=1.0/(exp((9.0-ee)/3.5)+1.0);
            //    Temp=hbarc*((pow(0.087*ee,0.25)*Aper  \
            //	  +pow(0.102*ee,0.25)*(1.0-Aper) \
            //	  +0.15*(1.0-Bper))*Cper		\
            //      +pow(0.15*ee,0.15)*(1-Cper));
            //    
            //    if(Temp<0.165)
            //      {}
            break;

        case EOSL2:
            Temp = eosl.T(eps);
            break;
        default:
            Temp=(eps+P(eps,bd))/S(eps,bd);
    }
    return Temp;
}

double CEos::S(CD &eps, CD &bd)
{
    double Sd; //entropy density
    double gt=169.0/4.0;         //Total freedom of Quarks and Gluon
    double ee;
    double a=(gt*pi*pi)/90.0;

    CD eosl_e1=0.1270769021427449;
    CD eosl_e2=0.446707952467404;
    CD eosl_e3=1.9402832534193788;
    CD eosl_e4=3.7292474570977285;
    switch(IEOS){
        ee=eps/hbarc;
        case EOSI:
        break;
        case EOSQ:
        break;
        case SMEOSQ:
        break;
        case EOSL:
        if(eps<eosl_e1)Sd=12.2304*pow(eps,1.16849);
        else if(eps>=eosl_e1 && eps<eosl_e2)Sd=11.9279*pow(eps,1.15635);
        else if(eps>=eosl_e2 && eps<eosl_e3)Sd=0.0580578+11.833*pow(eps,1.16187);
        else if(eps>=eosl_e3 && eps<eosl_e4)Sd=18.202*eps-62.021814-4.85479*exp(-2.72407*1.0E-11*pow(eps,4.54886))\
            +65.1272*pow(eps,-0.128012)*exp(-0.00369624*pow(eps,1.18735))-4.75253*pow(eps,-1.18423);
        else if(eps>eosl_e4)Sd=18.202*eps-62.021814-4.85479*exp(-2.72407*1.0E-11*pow(eps,4.54886))\
            +65.1272*pow(eps,-0.128012)*exp(-0.00369624*pow(eps,1.18735));
        break;
        case EOSL2:
        Sd=eosl.S(eps);
        break;
        default:
        Sd=0.0;
    }
    return Sd;
}

//#define TESTEOS
#ifdef TESTEOS
int main()
{
    const int N=300;
    CD de=600/double(N);
    ofstream fP("P_eosl.dat");
    ofstream fP_eosi("P_eosi.dat");
    ofstream fP_eosq("P_eosq.dat");
    ofstream fS("S_eosl.dat");
    ofstream fT("T_eosl.dat");
    ofstream fP_eosl2("P_eosl2.dat");	

    CEos eos(EOSL);           //Create a EOSL object
    CEos eosi(EOSI);          //Create a EOSI object
    CEos eosq(SMEOSQ);
    CEos eosl(EOSL2);

    for(int i=0; i<N; i++){
        double eps=i*de;
        double bd=0.0;
        fP<<eps<<' '<<eos.P(eps,bd)<<endl;
        fP_eosi<<eps<<' '<<eosi.P(eps,bd)<<endl;
        fP_eosq<<eps<<' '<<eosq.P(eps,bd)<<endl;
        fS<<eps<<' '<<eos.S(eps,bd)<<endl;
        fT<<eps<<' '<<eos.T(eps,bd)<<endl;

        fP_eosl2<<eps<<' '<<eosl.P(eps,bd)<<endl;
    }
    cout<<"e=0.12, T="<<eosl.T(0.120,0.0)<<endl;
    cout<<"e=0.11, T="<<eosl.T(0.110,0.0)<<endl;
    cout<<"e=0.075, T="<<eosl.T(0.075,0.0)<<endl;
    fP.close();
    fS.close();
    fT.close();
    fP_eosl2.close();
}
#endif
