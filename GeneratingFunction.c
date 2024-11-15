#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
const long long p=1000000007;
long long n,ans;
long long X[200000]={0,1},cnt[200000]={0,1,2,1},tmp[200000],ttt[200000];
void polymul(long long *a,long long *b,long long *c,long long n)
{
    long long *ans=(long long*)calloc(n+1,sizeof(long long));
    for(long long i=0;i<=n;++i)
        for(long long j=0;j<=i;++j)
            ans[i]=(ans[i]+a[j]*b[i-j]%p)%p;    //ans[i]=a[0]*b[i]+a[1]*b[i-1]+...+a[i]*b[0]
    for(long long i=0;i<=n;++i)
        c[i]=ans[i];
}
int main()
{
    scanf("%lld",&n);
    ans=(ans+cnt[n])%p;
    for(long long h=1;h<=log2(n+1);++h)
    {                               //cnt=T_h(x)
        polymul(X,cnt,tmp,n);       //tmp=xT_h(x)
        polymul(tmp,cnt,ttt,n);     //ttt=xT_h(x)T_h(x)
        tmp[0]=1;                   //tmp=xT_h(x)+1
        polymul(tmp,tmp,tmp,n);     //tmp=(xT_h(x)+1)^2 
        polymul(ttt,tmp,cnt,n);     //T_{h+1}(x)=xT_h(x)T_h(x)(xT_h(x)+1)^2
        ans=(ans+cnt[n])%p;
    }                               //cnt=T_{h+1}(x)
    printf("%lld\n",ans);
    return 0;
}