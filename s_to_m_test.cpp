#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void error(int mark)
{
    cout<<"Error: "<<mark<<endl;
}

vector<int> s_to_m(string s)
{
    vector<int> iter;
    int digit = 0;
    int num = 0;
    bool in_bracket = 0;
    bool in_bracket_out = 0;

    for(int i=0;i<s.length();i++)
    {
        if(s[i] == ' ')
            continue;
        else if(s[i]=='{' && !in_bracket_out && !in_bracket)
        {
            in_bracket_out = 1;
            continue;                
        }
        else if(s[i]=='}' && in_bracket_out && !in_bracket)
        {
            in_bracket_out = 0;
            continue;    
        }
        else if(s[i]==',' && !in_bracket && in_bracket_out)
        {
            continue;
        }

        if(in_bracket_out)
        {
            if(s[i]=='{' && !in_bracket)
            {               
                in_bracket = 1; 
                continue;
            } 
            else if(s[i]=='}' && in_bracket)
            {
                iter.push_back(num);
                digit = 0;
                num = 0;
                in_bracket = 0;
                continue;
            }
            else if(s[i]==',' && in_bracket)
            {
                iter.push_back(num);
                digit = 0;
                num = 0;
            }
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
            else
            {
                error(2);
                break;
            }
        } 
    }
    return iter;
}


int main()
{
    string s;
    s = "3 {{1,3,5,23} , {23,54,2,3}} ";

    vector<int> a = s_to_m(s);
    for(int i=0;i<a.size();i++)
    {
        cout<<a[i]<<" ";
    }
    cout<<endl;


    return 0;
}

