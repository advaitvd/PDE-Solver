#pragma once
using namespace std;
#include<iostream>
#include<unordered_map>
#include<stack>
#include "AD.h"

template<class T>
class charList{
	private:
		string expression;
		vector<T> nums;
		vector<AD<T>> indVars;
	public:
		charList();
		charList(string str);
		string infix_to_postfix(unordered_map<string, AD<T>>&vars);
		AD<T> evaluate(int varCount);
		void display();
		int precedence(char c);
};

template<class T>
charList<T>::charList():expression(""){}

template<class T>
charList<T>::charList(string str)
	:expression(str){}

template<class T>
int charList<T>::precedence(char c){
	if(c=='^'){
		return 5;
	}else if(c=='/'||c=='*'){
		return 4;
	}else if(c=='+'||c=='-'){
		return 3;
	}
}

template<class T>
string charList<T>::infix_to_postfix(unordered_map<string, AD<T>>&vars){
	stack<char> stk;
	stack<char> parentheses;
	vector<int> count;
	string res="";
	for(int i=0;i<expression.size();i++){
		if(isalpha(expression[i])||expression[i]=='_'){
			string temp="";
			while(i<expression.size()&&(isalnum(expression[i])||expression[i]=='_')){
				temp+=expression[i];
				i++;
			}
			i--;
			if(vars.find(temp)!=vars.end()){
				res+="#";
				indVars.push_back(vars[temp]);
			}else{
				res+=(temp+"(");
				count.push_back(parentheses.size()+1);
			}
		}else if(isdigit(expression[i])){
			res+='$';
			string temp="";
			while(i<expression.size() && (isdigit(expression[i])||expression[i]=='.')){
				temp+=expression[i];
				i++;
			}
			i--;
			nums.push_back(stod(temp));
		}else if(expression[i]=='('||stk.empty()||(stk.top()=='('&&expression[i]!=')')){
			stk.push(expression[i]);
			if(expression[i]=='('){
				parentheses.push('(');
			}
		}else if(expression[i]==')'){
			while(stk.top()!='('){
				res+=stk.top();
				stk.pop();
			}
			stk.pop();
			if(!count.empty()&&*(count.end()-1)==parentheses.size()){
				res+=')';
				count.pop_back();
			}
			parentheses.pop();
		}else{
			if(precedence(expression[i])>precedence(stk.top())){
				stk.push(expression[i]);
			}else{
				while(!stk.empty()&&stk.top()!='('&&precedence(stk.top())>=precedence(expression[i])){
					res+=stk.top();
					stk.pop();
				}
				stk.push(expression[i]);
			}
		}
	}
	while(!stk.empty()){
		res+=stk.top();
		stk.pop();
	}
	this->expression=res;
	return res;
}

template<class T>
AD<T> charList<T>::evaluate(int varCount){
	stack<string> stk_func;
	stack<AD<T>> stk_main;
	
	int k=0;
	int mm=0;
	for(int i=0;i<expression.size();i++){
		if(expression[i]=='('){
			continue;
		}else if(expression[i]=='$'){
			AD<T> tmp;
			tmp.set_f(nums[k++]);
			tmp.set_df_const(varCount);
			stk_main.push(tmp);
		}else if(expression[i]=='#'){
			stk_main.push(indVars[mm++]);
		}else if(isalpha(expression[i])){
			string fn="";
			while(i<expression.size()&&isalpha(expression[i])){
				fn+=expression[i];
				i++;
			}
			i--;
			stk_func.push(fn);
		}else if(expression[i]==')'){
			AD<T> arg=stk_main.top();
			stk_main.pop();
			if(stk_func.top()=="sin"){
				stk_main.push(sin(arg));
			}else if(stk_func.top()=="cos"){
				stk_main.push(cos(arg));
			}else if(stk_func.top()=="tan"){
				stk_main.push(tan(arg));
			}else if(stk_func.top()=="cosec"){
				stk_main.push(cosec(arg));
			}else if(stk_func.top()=="sec"){
				stk_main.push(sec(arg));
			}else if(stk_func.top()=="cot"){
				stk_main.push(cot(arg));
			}else if(stk_func.top()=="arccos"){
				stk_main.push(arccos(arg));
			}else if(stk_func.top()=="arcsin"){
				stk_main.push(arcsin(arg));
			}else if(stk_func.top()=="arctan"){
				stk_main.push(arctan(arg));
			}else if(stk_func.top()=="sinh"){
				stk_main.push(sinh(arg));
			}else if(stk_func.top()=="cosh"){
				stk_main.push(cosh(arg));
			}else if(stk_func.top()=="tanh"){
				stk_main.push(tanh(arg));
			}else if(stk_func.top()=="log"){
				stk_main.push(log(arg));
			}else if(stk_func.top()=="exp"){
				stk_main.push(exp(arg));
			}else if(stk_func.top()=="abs"){
				stk_main.push(abs(arg));
			}
			stk_func.pop();
		}else{
			if(expression[i]=='^'){
				AD<T> b=stk_main.top();
				stk_main.pop();
				AD<T> a=stk_main.top();
				stk_main.pop();
				stk_main.push(a^b);
			}else if(expression[i]=='/'){
				AD<T> b=stk_main.top();
				stk_main.pop();
				AD<T> a=stk_main.top();
				stk_main.pop();
				stk_main.push(a/b);
			}else if(expression[i]=='*'){
				AD<T> b=stk_main.top();
				stk_main.pop();
				AD<T> a=stk_main.top();
				stk_main.pop();
				stk_main.push(a*b);
			}else if(expression[i]=='+'){
				AD<T> b=stk_main.top();
				stk_main.pop();
				AD<T> a=stk_main.top();
				stk_main.pop();
				stk_main.push(a+b);
			}else if(expression[i]=='-'){
				AD<T> b=stk_main.top();
				stk_main.pop();
				AD<T> a=stk_main.top();
				stk_main.pop();
				stk_main.push(a-b);
			}
		}
	}
	
	return stk_main.top();
}

template<class T>
void charList<T>::display(){
	for(int i=0;i<expression.size();i++){
		cout<<expression[i];
	}cout<<endl;
}

