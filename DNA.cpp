#include <bits/stdc++.h>
#include <stdio.h>
#include <iostream>
#include <string>

using namespace std;

 
// function to find out the minimum penalty
void getMinimumPenalty(char x[], char y[], int pxy, int pgap)
{
    int i, j; // initialising variables
     
    //int m = x.length(); // length of gene1
    //int n = y.length(); // length of gene2
    int m = 80, n = 80;
     
   // table for storing optimal substructure answers
    int dp[n+m+1][n+m+1] = {0};
 
    // initialising the table
    for (i = 0; i <= (n+m); i++)
    {
        dp[i][0] = i * pgap;
        dp[0][i] = i * pgap;
    }   
 
    // calculating the minimum penalty
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            if (x[i - 1] == y[j - 1])
            {
                dp[i][j] = dp[i - 1][j - 1];
            }
            else
            {
                dp[i][j] = min({dp[i - 1][j - 1] + pxy ,
                                dp[i - 1][j] + pgap    ,
                                dp[i][j - 1] + pgap    });
            }
        }
    }
 
    // Reconstructing the solution
    int l = n + m; // maximum possible length
     
    i = m; j = n;
     
    int xpos = l;
    int ypos = l;
 
    // Final answers for the respective strings
    int xans[l+1], yans[l+1];
     
    while ( !(i == 0 || j == 0))
    {
        if (x[i - 1] == y[j - 1])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--; j--;
        }
        else if (dp[i - 1][j - 1] + pxy == dp[i][j])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--; j--;
        }
        else if (dp[i - 1][j] + pgap == dp[i][j])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)'_';
            i--;
        }
        else if (dp[i][j - 1] + pgap == dp[i][j])
        {
            xans[xpos--] = (int)'_';
            yans[ypos--] = (int)y[j - 1];
            j--;
        }
    }
    while (xpos > 0)
    {
        if (i > 0) xans[xpos--] = (int)x[--i];
        else xans[xpos--] = (int)'_';
    }
    while (ypos > 0)
    {
        if (j > 0) yans[ypos--] = (int)y[--j];
        else yans[ypos--] = (int)'_';
    }
 
    // Since we have assumed the answer to be n+m long,
    // we need to remove the extra gaps in the starting
    // id represents the index from which the arrays
    // xans, yans are useful
    int id = 1;
    for (i = l; i >= 1; i--)
    {
        if ((char)yans[i] == '_' && (char)xans[i] == '_')
        {
            id = i + 1;
            break;
        }
    }
 
    // Printing the final answer
    cout << "Minimum Penalty in aligning the genes = ";
    cout << dp[m][n] << "\n";
    cout << "The aligned genes are :\n";
    for (i = id; i <= l; i++)
    {
        cout<<(char)xans[i];
    }
    cout << "\n";
    for (i = id; i <= l; i++)
    {
        cout << (char)yans[i];
    }
    return;
}




int main(){
    // initialising penalties of different types
    int misMatchPenalty = 3;
    int i, gapPenalty = 2;
    int len, row, chk = 1;
    FILE *gene1, *gene2;
    gene1 = fopen("GCF_000009045.1.fna", "r");
    gene2 = fopen("GCF_000242855.2.fna", "r");
    char a[90], b[90];
    if(gene1 && gene2 == NULL)
    {
              exit(1);
    }
    else
    {
        while(chk == 1)
        {
            fscanf(gene1,"%s",a);
            fscanf(gene2,"%s",b);
            if(!feof(gene1) && !feof(gene2))
            {
                getMinimumPenalty(a, b, misMatchPenalty, gapPenalty);
            }
            else
            chk = 0;
        }
    }
    fclose(gene1);
    fclose(gene2);
    return 0;
}