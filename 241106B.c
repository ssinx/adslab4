#include<stdio.h>
int n;
long long ans;
const long long p=1e9+7;
long long red[1000][1000];      //color[i][j] is the number of R-B trees with i internal nodes and j black nodes, 
long long black[1000][1000];    //with the root being certain color
int main()
{
    /* Initialize */
    black[0][0]=1;
    
    /* Input */
    scanf("%d",&n);

    /* Compute */
    for(int i=1;i<=n;++i)
        for(int j=0;j<=i;++j)
            for(int k=0;k<i;++k)
            {
                if(j>0)
                    black[i][j]=(black[i][j]+(black[k][j-1]+red[k][j-1])*(black[i-k-1][j-1]+red[i-k-1][j-1])%p)%p;
                    // the son of a black root can be either black or red
                red[i][j]=(red[i][j]+black[k][j]*black[i-k-1][j]%p)%p;
                    // but the son of a red root must be black
            }
    
    /* Sum up */
    for(int i=1;i<=n;++i)
        ans=(ans+black[n][i])%p;
    
    /* Output */
    printf("%lld\n",ans);
    
    return 0;
}