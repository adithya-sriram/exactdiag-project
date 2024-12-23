#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "op.h"
#include "util.h"

using namespace std;


namespace pyED {

    SingleSiteOp::SingleSiteOp(int n, int op) {
        this->site = n;
        this->optype = op;
        if (op <= 8)
            this->fermionic = true;
        else
            this->fermionic = false;
    }
    string SingleSiteOp::to_string() {
        int n = (this->site);
        switch(this->optype) {
            case 1:
                return "Sx(" + std::to_string(n) + ")"; 
            case 2:
                return "Sy(" + std::to_string(n) + ")"; 
            case 3:
                return "Sz(" + std::to_string(n) + ")"; 
            case 9:
                return "b(" + std::to_string(n) + ")";
            case 10:
                return "bD(" + std::to_string(n) + ")";  
        }
    }
    
    Operator::Operator() { }
    Operator Operator::operator+(Operator& o1) {
        Operator newOp;
        newOp.ops = this->ops;
        for(int i = 0; i < o1.ops.size(); i++)
            newOp.ops.push_back(o1.ops[i]);
        return newOp;
    }
    Operator Operator::operator-(Operator& o1) {
        Operator newOp;
        newOp.ops = this->ops;
        for(int i = 0; i < o1.ops.size(); i++) {
            newOp.ops.push_back(o1.ops[i]);
            complex<double> amp = newOp.ops[newOp.ops.size() - 1].amplitude;
            newOp.ops[newOp.ops.size() - 1].amplitude.real(amp.real() * -1);
            newOp.ops[newOp.ops.size() - 1].amplitude.imag(amp.imag() * -1);
        }
        return newOp;
    }
    Operator Operator::operator*(Operator& o1) {
        Operator newOp;
        newOp.ops = this->ops;
        for (int i = 0; i < o1.ops.size(); i++) {
            for (int j = 0; j < newOp.ops.size(); j++) {
                for (int k = o1.ops[i].terms.size() - 1; k >= 0; k--) {
                    newOp.ops[j].terms.insert(newOp.ops[j].terms.begin(), o1.ops[i].terms[k]);
                }
                newOp.ops[j].amplitude *= o1.ops[i].amplitude;
            }
        
        }
        return newOp;
    }
    Operator Operator::operator*(complex<double> amp1) {
        Operator newOp;
        newOp.ops = this->ops;
        for (int i = 0; i < newOp.ops.size(); i++) {
            newOp.ops[i].amplitude *= amp1;
        }
        return newOp;
    }
    Operator& Operator::operator+=(Operator& o1) {
        for(int i = 0; i < o1.ops.size(); i++)
            this->ops.push_back(o1.ops[i]);
        return *this;
    }
    Operator& Operator::operator*=(Operator& o1) {
        for (int i = 0; i < o1.ops.size(); i++) {
            for (int j = 0; j < this->ops.size(); j++) {
                for (int k = o1.ops[i].terms.size() - 1; k >= 0; k--) {
                    this->ops[j].terms.insert(this->ops[j].terms.begin(), o1.ops[i].terms[k]);
                }
                this->ops[j].amplitude *= o1.ops[i].amplitude;
            }
        
        }
        return *this;
    }
    Operator& Operator::operator*=(complex<double> amp1) {
        for (int i = 0; i < this->ops.size(); i++) {
            this->ops[i].amplitude *= amp1;
        }
        return *this;
    }
    Operator Operator::power(int n) {
        Operator newOp;
        newOp.ops = this->ops;
        for (int i = 0; i < n; i++) {
            newOp *= *this;
        }
        return newOp;
    }
    // Currently broken
    Operator Operator::dag() {
        Operator newOp;
        for (int i = 0; i < this->ops.size(); i++){
            OperatorTerm newTerm;
            newTerm.terms = vector<SingleSiteOp> (this->ops[i].terms.rbegin(), this->ops[i].terms.rend());
            newOp.ops.push_back( newTerm );
            for (int j = 0; j < newOp.ops[i].terms.size(); j++) {
                switch(newOp.ops[i].terms[j].optype) {
                    case(1) : newOp.ops[i].terms[j].optype = 1;
                    case(2) : newOp.ops[i].terms[j].optype = 2;
                    case(3) : newOp.ops[i].terms[j].optype = 3;
                    case(4) : newOp.ops[i].terms[j].optype = 5;
                    case(5) : newOp.ops[i].terms[j].optype = 4;
                }
            }
            newOp.ops[i].amplitude = conj(newOp.ops[i].amplitude);
        }
        return newOp;
    }

    string Operator::to_string() {
        string obj = "";
        for (int i = 0; i < this->ops.size(); i++) {
            complex<double> amp = (this->ops[i]).amplitude;
            obj += "(" + prettyPrintComplex(amp) + ")";
            for (int j = 0; j < this->ops[i].terms.size(); j++) {
                obj += this->ops[i].terms[j].to_string();
                if (j != this->ops[i].terms.size() - 1)
                    obj += " * ";
            }
            if (i != this->ops.size() - 1)
                obj += " + ";
        }
        return obj;
    }


    Operator Sx(int n) {
        OperatorTerm X;
        X.terms.push_back( SingleSiteOp(n, 1) );
        X.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(X);
        return newOp;
    }
    Operator Sy(int n) {
        OperatorTerm Y;
        Y.terms.push_back( SingleSiteOp(n, 2) );
        Y.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(Y);
        return newOp;
    }
    Operator Sz(int n) {
        OperatorTerm Z;
        Z.terms.push_back( SingleSiteOp(n, 3) );
        Z.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(Z);
        return newOp;
    }
    Operator Sp(int n) {
        OperatorTerm Sp;
        Sp.terms.push_back( SingleSiteOp(n, 4) );
        Sp.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(Sp);
        return newOp;
    }
    Operator Sm(int n) {
        OperatorTerm Sm;
        Sm.terms.push_back( SingleSiteOp(n, 5) );
        Sm.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(Sm);
        return newOp;
    }
    Operator b(int n) {
        OperatorTerm b;
        b.terms.push_back( SingleSiteOp(n, 9) );
        b.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(b);
        return newOp;
    }

    Operator bD(int n) {
        OperatorTerm bD;
        bD.terms.push_back( SingleSiteOp(n, 10) );
        bD.amplitude.real(1);
        Operator newOp;
        newOp.ops.push_back(bD);
        return newOp;
    }


}