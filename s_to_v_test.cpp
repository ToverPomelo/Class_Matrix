#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void error(int mark)
{
    cout<<mark<<endl;
}

vector<int> s_to_v(string s)
{
    vector<int> iter;
    int digit = 0;
    int num = 0;
    bool in_bracket = 0;
//    if(s[0]!='{' || s[s.length()-1]!='}') 
//        error(1);
    for(int i=0;i<s.length();i++)
    {
        if(s[i] == ' ')
            continue;
        else if(s[i] == '{')
            in_bracket = 1; 
        else if((s[i]==',' || s[i]=='}') && in_bracket)
        {
            iter.push_back(num);
            digit = 0;
            num = 0;
        }
        else if(s[i] == '}')
            in_bracket = 0;
        else if(s[i]>='0' && s[i]<='9' && in_bracket)
        {
            int iter_num = (int)s[i] - (int)'0'; 
            num = num*pow(10,digit) + iter_num;
            digit++;
        }
        else if(in_bracket)
        {
            error(1);
            break;
        }
    }
    return iter;
}

int main()
{
    string s;
    s = "3 { 1, 3, 5, 23 } ";

    vector<int> a = s_to_v(s);
    for(int i=0;i<a.size();i++)
    {
        cout<<a[i]<<" ";
    }
    cout<<endl;


    return 0;
}

