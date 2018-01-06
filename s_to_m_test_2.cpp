#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void error(int mark)
{
    cout<<"Error: "<<mark<<endl;
}

vector<double> s_to_m(string s)
{
    vector<double> iter;
    int digit_left = 0;
    int digit_right = 0;
    double num = 0;
    bool in_bracket = 0;
    bool in_bracket_out = 0;
    bool is_double = 0;

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
                digit_left = 0;
                digit_right = 0;
                is_double = 0;
                num = 0;
                in_bracket = 0;
                continue;
            }
            else if(s[i]==',' && in_bracket)
            {
                iter.push_back(num);
                digit_left = 0;
                digit_right = 0;
                is_double = 0;
                num = 0;
            }
            else if(s[i]=='.' && in_bracket)
            {
                is_double = 1;
                continue;
            }
            else if(s[i]>='0' && s[i]<='9' && in_bracket && !is_double)
            {
                int iter_num = (int)s[i] - (int)'0'; 
                num = num*pow(10,digit_left) + iter_num;
                digit_left++;
            }
            else if(s[i]>='0' && s[i]<='9' && in_bracket && is_double)
            {
                digit_right++;
                int iter_num = (int)s[i] - (int)'0'; 
                num = num + iter_num*pow(10,(-1)*digit_right);
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
    s = "3 {{1,3,5.124,23.24} , {23,54.45,2,3}} ";

    vector<double> a = s_to_m(s);
    for(int i=0;i<a.size();i++)
    {
        cout<<a[i]<<" ";
    }
    cout<<endl;


    return 0;
}

