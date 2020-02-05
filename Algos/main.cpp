#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <windows.h>
#include <bits/stdc++.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>


using namespace std;

int Brute(string Table, string pattern, int n, int b)
{
    int match=0;
    for (int i=0; i<n; i++)
    {
        int found=0;
        for (int j=0; j<b; j++)
        {
            if(Table[i+j]==pattern[j])
            {
                found=1;
            }
            else
            {
                found=0;
                break;
            }
        }
        if(found==1)
           match++;
    }
    return match;
}

//*****************************************************************************************

int Sunday(string Table, string pattern, int n, int b)
{
    bool newTable[256] = {false};
    int match=0;
    for (int i=0; i<b; i++)
        newTable[pattern[i]]= true;
    for (int i=0; i<n; i++)
    {
        int found=0;
        for (int j=0; j<b; j++)
        {
            if(Table[i+j]==pattern[j])
            {
                found=1;
            }
            else
            {
         if(newTable[Table[i+b]]==false)
            {
            i=i+b;
            }
                found=0;
                break;
            }
        }
        if(found==1)
            match++;
    }
    return match;
}

//*****************************************************************************************

int State(string pat, int siz, int state, int NumChar)
{
    if (state < siz && NumChar == pat[state])
        return state+1;
    int State_i, i;

    for (State_i = state; State_i > 0; State_i--)
    {
        if (pat[State_i-1] == NumChar)
        {
            for (i = 0; i < State_i-1; i++)
                if (pat[i] != pat[state-State_i+1+i])
                    break;
            if (i == State_i-1)
                return State_i;
        }
    }
    return 0;
}


void Compute(string pat, int siz, int **TF)
{
    int state, NumChar;
    for (state = 0; state <= siz; ++state)
        for (NumChar = 0; NumChar < 255; ++NumChar)
            TF[state][NumChar] = State(pat, siz, state, NumChar);
}


int FSM(string text, string pat, int n, int b)
{
    int match = 0;
    int** TF = new int*[b+1];
    for(int i = 0; i < b+1; ++i)
        TF[i] = new int[255];

    Compute(pat, b, TF);
    int i, state=0;
    for (i = 0; i < n; i++)
    {
        state = TF[state][text[i]];
        if (state == b)
        {
            match++;
        }
    }
    return match;
}


//*****************************************************************************************

unsigned const int W = 997;
unsigned const int S = 128;


unsigned int hashing(string s, int len, int j=0)
{
    unsigned  int h =0;
    for (int k=j; k<len; k++)
    {
        h=((h*S)%W+s[k])%W;
    }
    return h;
}

unsigned int rehash ( unsigned int Spower, char oldc, char newc, unsigned int oldhash)
{
    oldhash = (((oldhash - (oldc * Spower)%W+W)%W)*S%W+newc)%W;
    return oldhash;
}

size_t Power(size_t n, size_t m)
{
    size_t res = 1;
    for (int i = 0 ; i < m; i++)
        res = (res * n )% W;
    return res;
}

int Hashy (string Fi, string Se, int sizFir, int sizSec)
{
    int found =0, i = 0, H1 =0, H2=0, match =0;
    size_t Spow = Power(S,sizSec-1);
    H2=hashing (Se, sizSec);
    for (i; i<=sizFir-sizSec; i++)
    {
        if (i==0)
        {
            H1=hashing(Fi, sizSec);
        }
        else
        {
                H1=rehash(Spow, Fi[i-1], Fi[i+sizSec-1], H1);
        }
        if(H1==H2)
        {
            match++;
            for (int j = 0; j < Se.size(); j++){
             if (Fi[i+j] != Se[j]){
                 match--;
                 break;
             }
        }
        }
    }
    return match;
}


//*****************************************************************************************



void PreProcKmp(string Pat, int n, int* lps)
{
    int i=1;
    int len =0;
    lps[0]=0;
    while (i < n)
    {
        if(Pat[i]==Pat[len])
        {
            len++;
            lps[i]=len;
            i++;
        }
        else
        {
            if(len!=0)
            {
                len=lps[len-1];
            }
            else
            {
                lps[i]=0;
                i++;
            }
        }
    }
}

int KMP (string Table, string pattern, int n, int b)
{
    int* lps = new int[b];
    int match=0;
PreProcKmp(pattern, b, lps);
int i=0, j=0;
while(i<n)
{
if(Table[i]==pattern[j])
{
    i++;j++;
}
if (j==b)
{
    match++;
    j=lps[j-1];
}
else
{
    if(i<n && Table[i]!=pattern[j])
    {
        if(j!=0)
        j=lps[j-1];
    else
        i=i+1;
    }
}
}
return match;
}

//******************************************

#include <iostream>
#include <string>

using namespace std;



void hashing2D(string *s , int *h, int hei, int len)
{
    for (int k=0; k<hei; k++)
    {
        for (int i=0; i<len; i++)
        h[k]=((h[k]*S)%W+s[k][i])%W;
    }

}

void rehash2D ( unsigned int Spower, string *s, int hei, int len, int *oldhash ,int pos)
{
    for(int i=0; i<hei; i++)
    {
    oldhash[i] = (((oldhash[i] -(s[i][pos] * Spower)%W+W)%W)*S%W+s[i][len+pos])%W;
    }
}

size_t PowerOk(size_t n, size_t m)
{
    size_t res = 1;
    for (int i = 0 ; i < m; i++)
        res = (res * n )% W;
    return res;
}

int SingleHash(int* Hash ,int len, int pos){
int H = 0;
for (int i = 0; i < len; i++){
    H = ((H*S) + Hash[i+pos])%W;
}
return H;
}

int SingleReHash(unsigned int Spower, int* Hash , int oldhash, int len, int pos){
  oldhash = (((oldhash -(Hash[pos] * Spower)%W+W)%W)*S%W+Hash[len+pos])%W;
  return oldhash;
}

int Hashy2d (string *Fi, string *Se, int sizFir, int sizSec)
{
    int found =0, i = 0, H1 =0, H2=0, match =0;
    int *HfT = new int[sizFir]();
    int *HfP = new int[sizSec]();
    size_t Spow = PowerOk(S,sizSec-1);
    hashing2D(Se , HfP, sizSec, sizSec);
    int PatHash = SingleHash(HfP,sizSec,0);
    int TextHash;
    for (int i = 0; i <= sizFir - sizSec; i++){

    if (i == 0)
        hashing2D(Fi , HfT, sizFir, sizSec);
    else
     rehash2D ( Spow, Fi, sizFir, sizSec, HfT , i-1);
    for (int j; j<=sizFir-sizSec; j++)
    {
        if (j == 0)
        TextHash = SingleHash(HfT, sizSec,0);
      else
        TextHash = SingleReHash(Spow,HfT,TextHash,sizSec,j - 1);
        if(TextHash==PatHash)
        {
            match++;
            for (int k = 0; k < sizSec; k++){
             if (Fi[j+k] != Se[k]){
                 match--;
                 break;
             }
        }
        }
        }
    }
    return match;
}





//*****************************************************************************************


int main()
{
//        ofstream myfile;
//ifstream fin("C:\\Users\\Nik_g\\OneDrive\\Documents\\study code\\Algos\\Text3.txt");
ofstream res("C:\\Users\\Nik_g\\OneDrive\\Documents\\study code\\Algos\\check.txt");
//string Table((istreambuf_iterator<char>(fin)), istreambuf_iterator<char>());


LARGE_INTEGER li;
QueryPerformanceFrequency(&li);
double PCFreq = li.QuadPart/1000.0;
int64_t startt;
int64_t endd;
std::int64_t taken1 =0;
std::int64_t taken2=0;
std::int64_t taken3=0;
std::int64_t taken4=0;
std::int64_t taken5=0;

int i=0;
    string Table = "qwe dasdfasczxvbsdbsdmasfmsaffsdfsdmgsdmfaskaskdasfjadghdsndsgjeqwejtqwtuydhfhafsj";
    string pattern = "qwe";
    int siz1 = Table.size();
    int siz2 = pattern.size();
  cout   << "Size of Table: "<< siz1 << "\t Size of pattern: " << siz2 << endl;
    if (siz2<=siz1)
    {
for (i; i<100; i++)
{
        QueryPerformanceCounter(&li);
startt = li.QuadPart;
    Brute(Table, pattern, siz1, siz2);
QueryPerformanceCounter(&li);
endd = li.QuadPart;
taken1 = endd - startt + taken1;


 QueryPerformanceCounter(&li);
startt = li.QuadPart;
  Hashy(Table, pattern, siz1, siz2);
QueryPerformanceCounter(&li);
endd = li.QuadPart;
taken2 = endd - startt + taken2;
QueryPerformanceCounter(&li);
startt = li.QuadPart;
  Sunday(Table, pattern, siz1, siz2);
QueryPerformanceCounter(&li);
endd = li.QuadPart;
taken3 = endd - startt + taken1;
QueryPerformanceCounter(&li);
startt = li.QuadPart;
   KMP(Table, pattern, siz1, siz2);
QueryPerformanceCounter(&li);
endd = li.QuadPart;
taken4 = endd - startt + taken1;

QueryPerformanceCounter(&li);
startt = li.QuadPart;
FSM(Table, pattern, siz1, siz2);
QueryPerformanceCounter(&li);
endd = li.QuadPart;
taken5 = endd - startt + taken1;

    }
    }
cout << FSM(Table, pattern, siz1, siz2) << endl;
cout << KMP(Table, pattern, siz1, siz2) << endl;
cout << Brute(Table, pattern, siz1, siz2) << endl;
cout << Hashy(Table, pattern, siz1, siz2) << endl;
cout << Sunday(Table, pattern, siz1, siz2) << endl;

cout <<" Brute: "<< std::setw(8) << std::fixed << " " << " " <<taken1/PCFreq/i<< std::endl;
cout <<" Carp: "<< std::setw(8) << std::fixed << " " << " " <<taken2/PCFreq/i<< std::endl;
cout <<" Sunday: "<< std::setw(8) << std::fixed << " " << " " <<taken3/PCFreq/i<< std::endl;
cout <<" KMP: "<< std::setw(8) << std::fixed << " " << " " <<taken4/PCFreq/i<< std::endl;
cout <<" FSM: "<< std::setw(8) << std::fixed << " " << " " <<taken5/PCFreq/i<< std::endl;








int TextLeng = 4;
int PatLeng = 2;
string text[TextLeng] = {"adoaa","aona","aaaa","aaaa"};
string pat[PatLeng] = {"do","on"};
cout<<Hashy2d(text,pat,TextLeng,PatLeng);







}

