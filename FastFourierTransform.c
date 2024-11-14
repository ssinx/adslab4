#include<stdio.h>
#include<math.h>
#include<string.h>
const double PI=3.14159265358979323846;
const long long p=1e9+7,M=1<<15,M2modp=(M*M)%p;//M2modp=M^2 mod p
long long n,ans,N;
double cnt[2000000]={0,1,2,1},tmp[2000000],ttt[2000000];
typedef struct{
    double Re,Im;
}Complex;
Complex omega(long long n,long long k,long long inv)
{
    Complex ans;
    ans.Re=cos(2*PI*k/n);
    if(inv)
        ans.Im=-sin(2*PI*k/n);
    else
        ans.Im=sin(2*PI*k/n);
    return ans;
}
Complex add(Complex a,Complex b)
{
    Complex ans;
    ans.Re=a.Re+b.Re;
    ans.Im=a.Im+b.Im;
    return ans;
}
Complex sub(Complex a,Complex b)
{
    Complex ans;
    ans.Re=a.Re-b.Re;
    ans.Im=a.Im-b.Im;
    return ans;
}
Complex mul(Complex a,Complex b)
{
    Complex ans;
    ans.Re=a.Re*b.Re-a.Im*b.Im;
    ans.Im=a.Re*b.Im+a.Im*b.Re;
    return ans;
}
Complex div(Complex a,Complex b)
{
    Complex ans;
    ans.Re=(a.Re*b.Re+a.Im*b.Im)/(b.Re*b.Re+b.Im*b.Im);
    ans.Im=(a.Im*b.Re-a.Re*b.Im)/(b.Re*b.Re+b.Im*b.Im);
    return ans;
}
void fft(Complex *a,long long n,long long inv)
{
    if(n==1)
        return;
    int f[n],L=log2(n);
    f[0]=0;
    for(int i=1;i<n;i++)            //calculate the bit reverse of i
        f[i]=(f[i>>1]>>1)|((i&1)<<(L-1));
    for(int i=0;i<n;i++)
        if(i<f[i])                  //swap a[i] and a[f[i]]
        {
            Complex t=a[i];
            a[i]=a[f[i]];
            a[f[i]]=t;
        }
    for(int mid=1;mid<n;mid<<=1)    //butterfly operation
    {
        Complex Wn=omega(n,n/(2*mid),inv);
        for(int R=mid<<1,j=0;j<n;j+=R)
        {
            Complex w;
            w.Re=1;
            w.Im=0;
            for(int k=0;k<mid;k++)
            {
                Complex x=a[j+k],y=mul(a[j+k+mid],w);
                a[j+k]=add(x,y);
                a[j+k+mid]=sub(x,y);
                w=mul(w,Wn);
            }
        }
    }
}
void polymul(double *a,double *b,double *c,long long n)
{
    Complex Cal[n],Cbl[n],Cah[n],Cbh[n];
    long long carrier;
    for(long long i=0;i<n;++i)
    {
        carrier=a[i];                   //a[i]=al[i]+M*ah[i]
        Cal[i].Re=carrier%M;
        Cal[i].Im=0;
        Cah[i].Re=carrier/M;
        Cah[i].Im=0;
        carrier=b[i];                   //b[i]=bl[i]+M*bh[i]
        Cbl[i].Re=carrier%M;
        Cbl[i].Im=0;
        Cbh[i].Re=carrier/M;
        Cbh[i].Im=0;
    }
    fft(Cal,n,0);                       //fft
    fft(Cbl,n,0);
    fft(Cah,n,0);
    fft(Cbh,n,0);
    Complex Ccll[n],Cchl[n],Cclh[n],Cchh[n];
    for(long long i=0;i<n;++i)          //a[i]*b[i]
    {                                   //=(al[i]+M*ah[i])*(bl[i]+M*bh[i])
        Ccll[i]=mul(Cal[i],Cbl[i]);     //=al[i]*bl[i]
        Cchl[i]=mul(Cah[i],Cbl[i]);     //+ah[i]*bl[i]*M
        Cclh[i]=mul(Cal[i],Cbh[i]);     //+al[i]*bh[i]*M
        Cchh[i]=mul(Cah[i],Cbh[i]);     //+ah[i]*bh[i]*M^2
    }
    fft(Ccll,n,1);                      //inverse fft
    fft(Cchl,n,1);
    fft(Cclh,n,1);
    fft(Cchh,n,1);
    long long ll,lh,hl,hh;
    for(long long i=0;i<n;++i)
    {
        ll=round(Ccll[i].Re/n);
        lh=round(Cclh[i].Re/n);
        hl=round(Cchl[i].Re/n);
        hh=round(Cchh[i].Re/n);
        c[i]=(ll%p+(hl+lh)%p*M%p+hh%p*M2modp%p)%p;
    }
}
void clear(double *a,long long n,long long N)   //clear a[n+1]~a[N]
{                                               //which are not used
    for(long long i=n;i<=N;++i)                 //to avoid bugs
        a[i]=0;
}
int main()
{
    scanf("%lld",&n);
    N=log2(n+1)+1;
    N=2<<N;                         //guarantee N>=2n
    long long carrier=cnt[n];
    ans=(ans+carrier)%p;
    for(long long h=1;h<=log2(n+1);++h)
    {                               //cnt=T_h(x)
        for(int i=N;i>=1;--i)       //tmp=xT_h(x)
            tmp[i]=cnt[i-1];
        tmp[0]=0;
        clear(tmp,n+1,N);
        polymul(tmp,cnt,ttt,N);     //ttt=xT_h(x)T_h(x)
        clear(ttt,n+1,N);
        tmp[0]=1;                   //tmp=xT_h(x)+1
        polymul(tmp,tmp,tmp,N);     //tmp=(xT_h(x)+1)^2 
        clear(tmp,n+1,N);
        polymul(ttt,tmp,cnt,N);     //T_{h+1}(x)=xT_h(x)T_h(x)(xT_h(x)+1)^2
        clear(cnt,n+1,N);
        carrier=cnt[n];
        ans=(ans+carrier)%p;
    }                               //cnt=T_{h+1}(x)
    printf("%lld\n",ans);
    return 0;
}