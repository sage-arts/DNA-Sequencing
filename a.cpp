#include<bits/stdc++.h>
using namespace std;

string get_hexamers_as_text(string word){
     int size=6,n=word.size();
     string hexamers_as_text="";
     for(int i=0;i<=n-size;i++){
          hexamers_as_text+=" "+word.substr(i,size);
     }
     return hexamers_as_text;
}

void fit(vector<string> human_texts,unordered_map<string,int>& word_idx){
     set<string> unique_words;
     string word;
     for(int i=0;i<human_texts.size();i++){
          stringstream ss(human_texts[i]);
          while(ss>>word){
               unique_words.insert(word);
          }
     }
     int i=0;
     for(auto it=unique_words.begin();it!=unique_words.end();it++){
          word_idx[*it]=i;
          i++;
     }
}

void transform(vector<string> human_texts,unordered_map<string,int> word_idx,int dim1,int dim2,vector<vector<int>>& _X){
     string word;
     for(int i=0;i<human_texts.size();i++){
          stringstream ss(human_texts[i]);
          while(ss>>word){
               _X[i][word_idx[word]]++;
          }
          i++;
     }
}

vector<double> get_mean(vector<vector<int>> mat){
     vector<double> mean(mat[0].size(),0);
     for(int i=0;i<mat.size();i++){
          for(int j=0;j<mean.size();j++){
               mean[j]+=mat[i][j];
          }
     }
     for(int i=0;i<mean.size();i++){
          mean[i]/=mat.size();
     }
     return mean;
}

int min_dist_classifier(vector<int> test,vector<vector<double>> u){
     int label;
     double min_dist=DBL_MAX,cur_dist;
     for(int i=0;i<u.size();i++){
          cur_dist=0;
          for(int j=0;j<test.size();j++){
               cur_dist+=pow(test[j]-u[i][j],2);
          }
          if(cur_dist<min_dist){
               min_dist=cur_dist;
               label=i;
          }
     }
     return label;
}

int main(){
     fstream file;
     string word,filename;
     filename="human_data.txt";
     file.open(filename.c_str());
     vector<string> human_texts;
     vector<int> y_data;
     int i=0;
     while(file>>word){
          if(i<2){
               i++;
               continue;
          } 
          if(i%2){
               y_data.push_back(stoi(word));
          }
          else{
               human_texts.push_back(get_hexamers_as_text(word));
          }
          i++;
     }
     unordered_map<string,int> word_idx;
     fit(human_texts,word_idx);
     int dim1=y_data.size(),dim2=word_idx.size();
     vector<vector<int>> _X(dim1,vector<int>(dim2,0));
     transform(human_texts,word_idx,dim1,dim2,_X);
     vector<vector<vector<int>>> X(7);
     for(int i=0;i<dim1;i++){
          X[y_data[i]].push_back(_X[i]);
     }
     vector<vector<int>> test_vectors;
     int samples=5;
     for(int i=0;i<X.size();i++){
          for(int j=0;j<samples;j++){
               test_vectors.push_back(X[i][j]);
          }
          X[i].erase(X[i].begin(),X[i].begin()+5);
     }
     vector<vector<double>> u;
     for(int i=0;i<X.size();i++){
          u.push_back(get_mean(X[i]));
     }
     for(int i=0;i<test_vectors.size();i++){
          int label=min_dist_classifier(test_vectors[i],u);
          cout<<label<<endl;
     }
     return 0;
}