//https://www.geeksforgeeks.org/top-20-dynamic-programming-interview-questions/
#include <bits/stdc++.h>
using namespace std;

// 1.Longest Common Subsequence

int lcs(string a, string b)
{
    if (a.length() == 0 || b.length() == 0)
        return 0;
    vector<vector<int>> dp(a.length() + 1, vector<int>(b.length() + 1, 0));
    //dp[i][j]=the maximum length subsequences formed by a[0]...a[i] and b[0]...b[i]
    for (int i = 1; i <= a.length(); i++)
    {
        for (int j = 1; j <= b.length(); j++)
        {
            if (a[i - 1] == b[j - 1])
            {
                //So dp[i][j] = the max length subsequence formed by a[0]...a[i-1] and b[0]...b[j]
                dp[i][j] = dp[i - 1][j - 1] + 1;
            }
            else
            {
                //dp[i][j]=Max of Max length subsequence formed by a[0]..a[i-1] and b[0]..b[j] or a[0]..a[i] and b[0]..b[j-1]
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    return dp[a.length()][b.length()];
}

// 2. Longest increasing Subseuence

int lis(vector<int> a)
{
    /*
    dp[i] represents longest increasing subseuence that ends with i
    */
    if (a.size() == 0)
        return 0;
    vector<int> dp(a.size(), 1);
    for (int i = 1; i < a.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (a[j] <= a[i])
            {
                dp[i] = max(dp[i], dp[j] + 1);
            }
        }
    }
    return *max_element(dp.begin(), dp.end());
}

// 3. Edit Distance

int editDistance(string a, string b)
{

    vector<vector<int>> dp(a.length() + 1, vector<int>(b.length() + 1, 0));

    for (int i = 0; i <= a.length(); i++)
    {
        for (int j = 0; j <= b.length(); j++)
        {
            if (i == 0)
            {
                //j insertions
                dp[i][j] = j;
            }
            else if (j == 0)
            {
                //i deletions
                dp[i][j] = i;
            }
            else if (a[i - 1] == b[j - 1])
            {
                dp[i][j] = dp[i - 1][j - 1];
            }
            else
            {
                //inserting  the character 1 + dp[i][j - 1]

                //deleting the character  1 + dp[i - 1][j]

                //replace the character 1 + dp[i - 1][j - 1]

                dp[i][j] = min(1 + dp[i][j - 1], min(1 + dp[i - 1][j], 1 + dp[i - 1][j - 1]));
            }
        }
    }
    return dp[a.length()][b.length()];
}

// 4.Minimum Partition

int minSub(vector<int> a)
{
    /*
        dp[i][j]=true if a[0]....a[i-1] sums to j
    */
    int sum = 0;
    int n = a.size();
    for (auto x : a)
        sum += x;
    vector<vector<bool>> dp(n + 1, vector<bool>(sum + 1, false));
    for (int i = 0; i <= n; i++)
    {
        dp[i][0] = true;
    }
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= sum; j++)
        {
            //Number is excluded
            dp[i][j] = dp[i - 1][j];
            //If Number is less than sum and included
            if (a[i - 1] <= j)
            {
                dp[i][j] = dp[i][j] || dp[i][j - a[i - 1]];
            }
        }
    }
    int diff = INT_MAX;
    for (int j = sum / 2; j >= 0; j--)
    {
        if (dp[n][j])
        {
            diff = min(diff, sum - 2 * j);
            break;
        }
    }
    return diff;
}

// 5 Ways to cover a distance

int countWays(int n)
{
    //dp[i] represents total number of ways to reach i

    int dp[n + 1];
    dp[0] = 1;
    if (n >= 1)
        dp[1] = 1;
    if (n >= 2)
        dp[2] = 2;
    for (int i = 3; i <= n; i++)
    {
        dp[i] = dp[i - 1] + dp[i - 2] + dp[i - 3];
    }
    return dp[n];
}

// 6 Longest Path in matrix

int recursive_longestPath(vector<vector<int>> &dp, vector<vector<int>> &a, int i, int j)
{
    if (!(i < dp.size() && i >= 0 && j >= 0 && j < dp[0].size()))
        return 0;

    if (dp[i][j] != -1)
        return dp[i][j];
    int r = 0, d = 0, l = 0, u = 0;
    if (j + 1 < dp[i].size() && a[i][j] + 1 == a[i][j + 1])
    {
        r = recursive_longestPath(dp, a, i, j + 1);
    }
    if (i + 1 < dp.size() && a[i][j] + 1 == a[i + 1][j])
    {
        d = recursive_longestPath(dp, a, i + 1, j);
    }
    if ((j - 1) >= 0 && a[i][j] + 1 == a[i][j - 1])
    {
        l = recursive_longestPath(dp, a, i, j - 1);
    }
    if ((i - 1) >= 0 && a[i][j] + 1 == a[i - 1][j])
    {
        u = recursive_longestPath(dp, a, i - 1, j);
    }
    dp[i][j] = 1 + max(r, max(d, max(l, u)));
    return dp[i][j];
}

int findLongestPath(vector<vector<int>> &a)
{
    /*
        dp[i][j] stores the longest path from i,j
        answer  is max of all dp[i][j]
    */

    if (a.size() == 0)
        return 0;
    vector<vector<int>> dp(a.size(), vector<int>(a[0].size(), -1));
    int ans = INT_MIN;
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < a[0].size(); j++)
        {
            if (dp[i][j] == -1)
                dp[i][j] = recursive_longestPath(dp, a, i, j);
            ans = max(ans, dp[i][j]);
        }
    }
    return ans;
}

// 7.Subset Sum

bool subsetSum(vector<int> a, int k)
{
    /*

    dp[i][j] represents whether a subset from A[0]....A[i] with sum j is possible or not

    */
    vector<vector<bool>> dp(a.size() + 1, vector<bool>(k + 1, false));

    for (int i = 0; i <= a.size(); i++)
    {
        dp[i][0] = true;
    }

    for (int i = 0; i <= k; i++)
    {
        dp[0][i] = false;
    }

    for (int i = 1; i <= a.size(); i++)
    {
        for (int j = 1; j <= k; j++)
        {
            if (a[i - 1] > j)
            {
                dp[i][j] = dp[i - 1][j];
            }
            else
            {
                dp[i][j] = dp[i - 1][j] || dp[i - 1][j - a[i - 1]];
            }
        }
    }

    return dp[a.size()][k];
}

// 8 Optimal Strategy for a game

int optimisegame(vector<int> a)
{
    // dp[i][j] represents the maximum score possible for a[i]....a[j]
    vector<vector<int>> dp(a.size(), vector<int>(a.size()));
    for (int l = 0; l < dp.size(); l++)
    {
        for (int i = 0, j = l; j < dp.size(); i++, j++)
        {
            int x = (i + 2 <= j) ? dp[i + 2][j] : 0;
            int y = (i + 1 <= (j - 1) && j - 1 >= 0) ? dp[i + 1][j - 1] : 0;
            int z = (j - 2 >= 0) ? dp[i][j - 2] : 0;
            dp[i][j] = max(a[i] + min(x, y), a[j] + min(y, z));
        }
    }
    return dp[0][dp.size() - 1];
}

// 9 0-1 Knapsack

int knapsack(vector<int> &v, vector<int> &w, int k)
{
    // dp[i][j] represents max value possible by taking a[0]..a[i] that weighs j

    vector<vector<int>> dp(v.size() + 1, vector<int>(k + 1, 0));
    for (int i = 1; i <= v.size(); i++)
    {
        for (int j = 1; j <= k; j++)
        {
            if (w[i - 1] <= j)
            {
                dp[i][j] = max(dp[i - 1][j], v[i - 1] + dp[i - 1][j - w[i - 1]]);
            }
            else
            {
                dp[i][j] = dp[i - 1][j];
            }
        }
    }
    return dp[v.size()][k];
}

// 11 Shortest Common Subsequence

int scs(string a, string b)
{
    // a.length()+b.length()-length of the longest common subsequence

    /*
        dp[i][j] represents length of the shortest common subsequence
        of a[0]...a[i] and b[0]...b[j]
    */
    vector<vector<int>> dp(a.length() + 1, vector<int>(b.length() + 1));

    for (int i = 0; i < dp.size(); i++)
    {
        for (int j = 0; j < dp[i].size(); j++)
        {
            if (i == 0)
            {
                dp[i][j] = j;
            }
            else if (j == 0)
            {
                dp[i][j] = i;
            }
            else if (a[i - 1] == b[j - 1])
            {
                dp[i][j] = 1 + dp[i - 1][j - 1];
            }
            else
            {
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    return dp[a.length()][b.length()];
}

// 12 Matrix Chain Multiplication

int matrixMultiplication(vector<int> p)
{
    /*
        dp[i][j] represents the min cost of multiplying matices A[i]..A[j]
    */

    vector<vector<int>> dp(p.size(), vector<int>(p.size(), INT_MAX));

    for (int i = 0; i < p.size(); i++)
    {
        dp[i][i] = 0;
    }

    for (int l = 2; l < p.size(); l++)
    {

        for (int i = 1; i <= (p.size() - l); i++)
        {
            int j = i + l - 1;
            for (int k = i; k < j; k++)
                dp[i][j] = min(dp[i][j], (dp[i][k] + dp[k + 1][j] + p[i - 1] * p[k] * p[j]));
        }
    }
    return dp[1][p.size() - 1];
}

// 13 Partition problem

bool partition(vector<int> a)
{

    int sum = 0;
    for (int i = 0; i < a.size(); i++)
    {
        sum += a[i];
    }
    // If sum is odd not possible
    if (sum % 2 != 0)
        return false;
    /*
        dp[i][j] represents whether the sum j is possible with a subset of a[0]...a[i]
    */
    vector<vector<bool>> dp(a.size() + 1, vector<bool>(sum));
    for (int i = 0; i <= a.size(); i++)
    {
        for (int j = 0; j <= sum; j++)
        {
            if (j == 0)
            {
                dp[i][j] = true;
            }
            else if (i == 0)
            {
                dp[i][j] = false;
            }
            else if (a[i - 1] <= j)
            {
                dp[i][j] = dp[i - 1][j] | dp[i - 1][j - a[i - 1]];
            }
            else
            {
                dp[i][j] = dp[i - 1][j];
            }
        }
    }
    return dp[a.size()][sum / 2];
}

// 14 Rod Cutting

int rodCutting(vector<int> p)
{
    // dp[i] represents max cost by cutting the rod of length i by possible costs
    vector<int> dp(p.size() + 1, 0);
    for (int l = 0; l <= p.size(); l++)
    {
        for (int j = 0; j < l; j++)
        {
            dp[l] = max(dp[l], p[j] + dp[l - j - 1]);
        }
    }
    return dp[p.size()];
}

// 15 Coin Change

int coinChange(int n, vector<int> a)
{

    /*
        dp[i][j] represents no of ways of getting sum j with i elements 
    */

    vector<vector<int>> dp(a.size() + 1, vector<int>(n + 1, 0));

    for (int i = 0; i <= a.size(); i++)
    {
        dp[i][0] = 1;
    }

    for (int i = 1; i < dp.size(); i++)
    {
        for (int j = 1; j < dp[0].size(); j++)
        {
            dp[i][j] = dp[i - 1][j];
            if ((j - a[i - 1]) >= 0)
                dp[i][j] += dp[i][j - a[i - 1]];
        }
    }

    return dp[a.size()][n];
}

// 16 Word break problem

bool wordBreak(string s, vector<string> d)
{
    /*
        dp[i] represents weather prefix s[0]..s[i] is there in dictionary or not
    */

    vector<bool> dp(s.length(), false);

    set<string> st;
    for (int i = 0; i < d.size(); i++)
    {
        st.insert(d[i]);
    }

    for (int l = 0; l < s.length(); l++)
    {
        string temp = s.substr(0, l + 1);
        if (!dp[l] && st.find(temp) != st.end())
            dp[l] = true;

        if (dp[l])
        {
            for (int j = 1; j < (s.length() - l); j++)
            {
                string temp = s.substr(l + 1, j);

                if (!dp[l + j] && st.find(temp) != st.end())
                {
                    dp[l + j] = true;
                }
                if ((l + j) == s.length() - 1 && dp[l + j])
                    return true;
            }
        }
    }
    return false;
}

// 17 Maximum Product When Cutting rope

int maxProduct(int n)
{
    if (n <= 1)
        return 0;

    vector<int> dp(n + 1, 0);
    dp[2] = 1;

    for (int l = 3; l <= n; l++)
    {
        for (int j = 1; j < (l / 2) + 1; j++)
        {
            dp[l] = max(dp[l], j * (l - j));
            if ((l - j) >= 2)
            {
                dp[l] = max(dp[l], j * dp[l - j]);
            }
        }
    }
    return dp[n];
}

// 18 Dice throw

int diceThrow(int n, int m, int x)
{

    /*
        dp[i][j] represents the sum obtained by throwing i dices
    */

    vector<vector<int>> dp(n + 1, vector<int>(x + 1, 0));

    for (int i = 1; i <= m && i <= x; i++)
    {
        dp[1][i] = 1;
    }

    for (int i = 2; i <= n; i++)
    {
        for (int j = i; j <= x; j++)
        {
            for (int k = 1; k <= m && k <= x; k++)
            {
                if ((j - k) >= 0 && dp[i - 1][j - k] != 0)
                {
                    dp[i][j] += dp[i - 1][j - k];
                }
            }
        }
    }
    return dp[n][x];
}

// 19 Box stacking problem
struct box
{
    int l, b, h;
    box(int p, int q, int r)
    {
        l = p;
        b = q;
        h = r;
    }
};

bool cmp(box p, box q)
{
    return ((p.l * p.b) > (q.l * q.b));
}

int maxHeight(vector<int> l, vector<int> b, vector<int> h)
{

    /*
        dp[i] represents max height possible with ith box on the top of the stack
    */

    vector<box> a;

    for (int i = 0; i < l.size(); i++)
    {
        a.push_back(box(max(l[i], b[i]), min(l[i], b[i]), h[i]));
        a.push_back(box(max(b[i], h[i]), min(b[i], h[i]), l[i]));
        a.push_back(box(max(l[i], h[i]), min(l[i], h[i]), b[i]));
    }
    sort(a.begin(), a.end(), cmp);

    vector<int> dp(a.size(), 1);

    for (int i = 0; i < a.size(); i++)
    {
        dp[i] = a[i].h;
    }

    for (int i = 1; i < a.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            if ((a[i].l < a[j].l) && (a[i].b < a[j].b) && (dp[i] < dp[j] + a[i].h))
            {
                dp[i] = dp[j] + a[i].h;
            }
        }
    }
    return *max_element(dp.begin(), dp.end());
}

// 20 Egg Breaking Problem

int eggBreak(int n, int f)
{
    /*
        dp[i][j] represents min no of trials required for i eggs and j floors

        Suppose  a egg is dropped from the k th floor (1<=k<=j) 
        then there are two ways egg breaks or doesn't breaks

        1.Egg breaks :
          so we have to check with (i-1) eggs in lower k-1 floors
        2.Egg doesn't break :
            so we have to check with eggs in the upper j-k floors
            
        dp[i][j] is the min among 1 + max(dp[i-1][k-1],d[i][j-k])
    */
    vector<vector<int>> dp(n + 1, vector<int>(f + 1, INT_MAX));

    for (int i = 0; i <= f; i++)
    {
        // 1 egg i floors i trials
        dp[1][i] = i;
        // 0 eggs i floors 0 trials
        dp[0][i] = 0;
    }
    for (int i = 0; i <= n; i++)
    {
        // i eggs 0 floors 0 trials
        dp[i][0] = 0;
        //i eggs 1 floor 1 trial
        dp[i][1] = 1;
    }
    for (int i = 2; i <= n; i++)
    {
        for (int j = 2; j <= f; j++)
        {

            for (int k = 1; k <= j; k++)
            {
                dp[i][j] = min(dp[i][j], 1 + max(dp[i - 1][k - 1], dp[i][j - k]));
            }
        }
    }
    return dp[n][f];
}

int main()
{
    cout << ":)" << endl;
    return 0;
}