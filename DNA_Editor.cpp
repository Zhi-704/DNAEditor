// CID: 01705802
// g++ a2.cpp -std=c++11
// gnFOXF1.fna
// https://www.youtube.com/watch?v=gH476CxJxfg, a song to work through the day

#include <fstream>
#include <iostream>
#include <ostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include "tictoc.h" //USE CTRL+F TO LOOK FOR "TIME", DID NOT COMPILE WHEN USED SO HAVE COMMENTED OUT

using namespace std;

string checker(string filebeingchecked);

// used for input validation of file
string checker(string user_in){
  bool inputvalid = true;

  cin.ignore();
  while (inputvalid){
    ifstream check(user_in, ios::in);
    if (check && (user_in.find(".fna")!=-1)) {
      inputvalid = false;
      cout  << "\nFound "
            << user_in
            << endl;
      check.close();
    }
    else if (!check || (user_in.find(".fna")==-1)) {
      cout  << "\nUnable to open file: "
            << user_in
            << "\nPlease specify a valid file name with the correct format\n"
            << "> ";
      getline(cin,user_in);
    }
    else {
      cout << "INPUT ERROR IN CHECKER FUNCTION";
    }
  }
  return user_in;
}

struct base{
  char nuc;
  base* forward = nullptr;
  base* backward = nullptr;
};

class DNAseq{
  private:
    string file_name;
    base* b_start = nullptr;
    base* b_end = nullptr;
    DNAseq* nextseq;
    int length = 0;
  public:
    DNAseq();
    DNAseq(string filename, DNAseq* old_object);
    void baseseq(char c);
    void getdata(string filename);
    DNAseq* getnextseq();
    bool validify(string basesbeingchecked);
    bool validify(int numberofbases);
    bool validify2(int numberofbases);
    int validify3(int numberofbases, int basepos);
    void search(string base, int choicemaker);
    void add(string base, int basepos, int choicemaker);
    void deletion(int basepos, int del_length);
    void replace(string base, int basepos, int del_length, int choicemaker);
    void printprevbases(base* currentbase);
    void printnextbases(base* currentbase);
    void save();
    void submenu();
};

// initialise DNAsequence to be nothing
DNAseq::DNAseq(){
  nextseq = nullptr;
}

// builds base sequence where each base points to the adjacent bases
void DNAseq::baseseq(char c){
  base* builder;
  if (b_start == nullptr && b_end == nullptr){
    b_start = new base;
    b_start->nuc = c;
    b_end = b_start;
  }
  else {
    builder = b_end;
    b_end = new base;
    b_end->nuc = c;
    b_end->backward = builder;
    builder->forward = b_end;
  }
}

// reads file and builds base sequence
void DNAseq::getdata(string file){
  ifstream fna;
  char c;
  string descriptor;

  fna.open(file,ios::in);
  getline(fna,descriptor);
  cout  << endl
        << "\nLoading file "
        << file
        << "...\n";
  while(!fna.eof()){
    fna.get(c);
    if (fna.peek()=='\n'){
      fna.ignore();
    }
    baseseq(c);
    length++;
  }
  fna.close();
  cout  << endl;
}

// DNA sequence points to old sequence or null if empty
DNAseq::DNAseq(string filename, DNAseq* old){
  file_name=filename;
  nextseq = old;
  getdata(filename);
}

// checks for the right bases being inputted
bool DNAseq::validify(string input){
  bool inputvalid = true;
  const string nuc = "ACGTNUKSYMWRBDHV-";
  istringstream check(input);
  char c;

  while (!check.eof()) {
    check.get(c);
    if (nuc.find(c) == -1) {
      inputvalid = false;
      return 0;
    }
  }
  if (inputvalid)
    return 1;
  else {
    cout << "\nVALIDIFY1b ERROR\n";
    return 0;
  }
}

// checks if right position was inputted
bool DNAseq::validify(int input){
  bool inputvalid = true;

  if (input <= length && input > 0) {
    inputvalid = false;
    return 0;
  }
  if (inputvalid) {
    return 1;
  }
  else {
    cout << "\nVALIDIFY1a ERROR\n";
    return 0;
  }
}

// checks for right value inputted for delete funciton
int DNAseq::validify3(int input, int basepos){
  bool inputvalid = true;
  int remaining = length-basepos;
  string baseinput;

  while (inputvalid) {
    if (basepos == length)  {
      if (input <= remaining+1 && input > 0) {
        inputvalid = false;
        break;
      }
      else if (input > remaining || input <= 0) {
      cout  << "\nPlease input a valid base pair length\n"
            << "Must be higher than 0 and lower than "
            << remaining+2
            << endl
            << "> ";
      getline(cin, baseinput);
      input = stoi(baseinput);
      }
      else {
        cout << "\nVALIDIFY3a ERROR\n";
        return 1;
      }
    }
    else {
      if (input <= remaining+1 && input > 0) {
        inputvalid = false;
        break;
      }
      else if (input >= remaining+2 || input <= 0) {
      cout  << "\nPlease input a valid base pair length\n"
            << "Must be higher than 0 and lower than "
            << remaining+2
            << endl
            << "> ";
      getline(cin, baseinput);
      input = stoi(baseinput);
      }
      else {
        cout << "\nVALIDIFY3b ERROR\n";
        return 1;
      }
    }
  }
  return input;
}

// checks if right position was inputted for add function
bool DNAseq::validify2(int input){
  bool inputvalid = true;

  if (input <= length+1 && input > 0) {
    inputvalid = false;
    return 0;
  }
  if (inputvalid) {
    return 1;
  }
  else {
    cout << "\nVALIDIFY2 ERROR\n";
    return 0;
  }
}

// retrieves the pointer to the next sequence in database
DNAseq* DNAseq::getnextseq(){
  return nextseq;
}

// prints the preceding 10 bases or stops if reaches null
void DNAseq::printprevbases(base* currentbase){
  base* prevbases;
  int sizecounter = 0;
  prevbases = currentbase;

  while ((prevbases->backward != nullptr) && sizecounter < 10) {
    prevbases = prevbases->backward;
    sizecounter++;
  }
  while (prevbases != currentbase){
    cout << prevbases->nuc;
    prevbases = prevbases->forward;
  }
}

// prints the next 10 bases or stops if reaches null
void DNAseq::printnextbases(base* currentbase){
  base* nextbases;
  int sizecounter = 0;
  nextbases = currentbase;

  while ((nextbases->forward != nullptr) && sizecounter < 10) {
    nextbases = nextbases->forward;
    cout << nextbases->nuc;
    sizecounter++;
  }
}

// searches for matching bases inside sequence with string input or fna file
void DNAseq::search(string baseinput, int choicemaker){
  base* finder1;
  base* finder2;
  int counter = 1;
  int sizecounter = 0;
  int match = 0;
  bool existence = false;

  if (choicemaker == 1){
    finder1 = b_start;
    finder2 = finder1->forward;
    sizecounter = baseinput.length();

    while (finder1 != nullptr) {
      if (finder1->nuc==baseinput[0] && baseinput.length()==1){// for psychopaths who write in one base
        cout  << endl
              << "Match #"
              << match
              << endl
              << "Base pair positions:\t["
              << counter
              << ":"
              << counter
              << "]\n"
              << "Base pair length:\t1\n"
              << "Prev 10 base pairs:\t";
        printprevbases(finder1);
        cout  << endl
              << "region of interest:\t"
              << baseinput
              << endl
              << "Next 10 base pairs:\t";
        printnextbases(finder1);
        cout  << endl;
        match++;
        existence = true;
        finder1 = finder1->forward;
        counter++;
      }
      else if (finder1->nuc==baseinput[0] && baseinput.length()>1) {
        finder2 = finder1->forward;
        for (int i = 1; i < baseinput.length(); i++) {
          if (finder2->nuc==baseinput[i]) {
            if (i == baseinput.length()-1) {
              cout  << endl
                    << "Match #"
                    << match
                    << endl
                    << "Base pair positions:\t["
                    << counter
                    << ":"
                    << counter+sizecounter-1
                    << "]\n"
                    << "Base pair length:\t"
                    << sizecounter
                    << endl
                    << "Prev 10 base pairs:\t";
              printprevbases(finder1);
              cout  << endl
                    << "region of interest:\t"
                    << baseinput
                    << endl
                    << "Next 10 base pairs:\t";
              printnextbases(finder2);
              cout  << endl;
              match++;
              existence = true;
            }
            finder2 = finder2->forward;
          }
          else
            break;
        }
        finder1 = finder1->forward;
        counter++;
      }
      else {
        finder1 = finder1->forward;
        counter++;
      }
    }
    if (!existence) {
      cout  << endl
            << "No matches were found\n";
    }
  }
  else if (choicemaker == 2 || choicemaker == 3) {
    DNAseq* tempseq;
    tempseq = new DNAseq(baseinput, tempseq);
    finder1 = b_start;
    finder2 = finder1->forward;
    sizecounter = tempseq->length;
    bool existence2 = false;
    base* printer = tempseq->b_start;

    while (finder1 != nullptr) {
      if (finder1->nuc==tempseq->b_start->nuc && sizecounter==1){// for psychopaths who write in one base IN A FNA FILE
        if (choicemaker == 3) {
          cout  << endl
                << "Found in DNA file:\t"
                << file_name;
        }
        cout  << endl
              << "Match #"
              << match
              << endl
              << "Base pair positions:\t["
              << counter
              << ":"
              << counter
              << "]\n"
              << "Base pair length:\t1\n"
              << "Prev 10 base pairs:\t";
        printprevbases(finder1);
        cout  << endl
              << "region of interest:\t";
        while (printer != nullptr){
          cout  << printer->nuc;
          printer = printer->forward;
        }
        cout  << endl
              << "Next 10 base pairs:\t";
        printnextbases(finder1);
        cout  << endl;
        match++;
        existence2 = true;
        finder1 = finder1->forward;
        counter++;
        printer = tempseq->b_start;
      }
      else if (finder1->nuc==tempseq->b_start->nuc && sizecounter>1) {
        finder2 = finder1->forward;
        base* looker;
        looker = tempseq->b_start->forward;
        for (int i = 1; i < sizecounter; i++) {
          if (finder2->nuc==looker->nuc) {
            if (i == sizecounter-1) {
              if (choicemaker == 3) {
                cout  << endl
                      << "Found in DNA file:\t"
                      << file_name;
              }
              cout  << endl
                    << "Match #"
                    << match
                    << endl
                    << "Base pair positions:\t["
                    << counter
                    << ":"
                    << counter+sizecounter-1
                    << "]\n"
                    << "Base pair length:\t"
                    << sizecounter
                    << endl
                    << "Prev 10 base pairs:\t";
              printprevbases(finder1);
              cout  << endl
                    << "region of interest:\t";
              while (printer != nullptr){
                cout  << printer->nuc;
                printer = printer->forward;
              }
              cout  << endl
                    << "Next 10 base pairs:\t";
              printnextbases(finder2);
              cout  << endl;
              match++;
              printer = tempseq->b_start;
              existence2 = true;
            }
            finder2 = finder2->forward;
            looker = looker->forward;
          }
          else
            break;
        }
        finder1 = finder1->forward;
        counter++;
      }
      else {
        finder1 = finder1->forward;
        counter++;
      }
    }
    if (!existence2) {
      cout  << endl
            << "No matches were found";
      if (choicemaker ==  3) {
        cout  << " in DNA file"
              << file_name;
      }
      cout << endl;
    }
  }
  else
   cout << "ERROR IN CHOICEMAKER SEARCH FUNCTION";
}

// add bases to specified base position with string input or fna file
void DNAseq::add(string baseinput, int basepos, int choicemaker) {
  base* finder1;
  base* finder2;
  base* printer;
  int sizecounter = 0;
  DNAseq* thead;

  if (choicemaker == 1) {
    char c;
    istringstream tempstream(baseinput);
    thead = new DNAseq();

    while(!tempstream.eof()){
      tempstream.get(c);
      if (tempstream.peek()=='\n'){
        tempstream.ignore();
      }
      thead->baseseq(c);
      thead->length++;
    }
  }
  else if (choicemaker == 2) {
    thead = new DNAseq(baseinput, thead);
  }
  else {
    cout << "\nCHOICEMAKER ERROR IN ADD FUNCTION\n";
  }
  finder1 = b_start;
  finder2 = finder1->forward;
  sizecounter = thead->length;
  printer = thead->b_start;
  if (basepos >= 2 && basepos <= length) {
    for (int i = 2; i < basepos; i++) {
      finder1 = finder1->forward;
    }
    finder2 = finder1->forward;
    thead->b_start->backward = finder1;
    finder1->forward = thead->b_start;
    thead->b_end->forward = finder2;
    finder2->backward = thead->b_end;
  }
  else if (basepos == 1) {
    finder2 = finder1;
    thead->b_end->forward = finder1;
    finder1->backward = thead->b_end;
  }
  else if (basepos == length+1) {
    for (int i = 2; i < basepos; i++) {
      finder1 = finder1->forward;
    }
    finder1->forward = thead->b_start;
    thead->b_start->backward = finder1;
    finder2 = thead->b_end->forward;
  }
  else {
    cout << "\nBASEPOS ERROR IN ADD FUNCTION\n";
  }
  cout  << endl
        << "Base pair positions:\t["
        << basepos
        << ":"
        << basepos+sizecounter-1
        << "]\n"
        << "Base pair length:\t"
        << sizecounter
        << endl
        << "Prev 10 base pairs:\t";
  printprevbases(thead->b_start);
  cout  << endl
        << "region of interest:\t";
  while (printer != finder2){
    cout  << printer->nuc;
    printer = printer->forward;
  }
  cout  << endl
        << "Next 10 base pairs:\t";
  printnextbases(thead->b_end);
  cout  << endl;
  length = length+sizecounter;
}

// delete a specified number of bases at specified location
void DNAseq::deletion(int basepos, int del_length){
  base* finder1;
  base* finder2;
  base* pdeleter;
  base* printer;
  finder1 = b_start;
  finder2 = b_start;

  for (int i = 1; i < basepos; i++) {
    finder1 = finder1->forward;
  }
  finder2 = finder1;
  for (int i = 0; i < del_length; i++) {
    finder2 = finder2->forward;
    if (finder2 == nullptr) {
      break;
    }
  }
  printer = finder1;
  cout  << endl
        << "DNA deletion information\n"
        << "Base pair positions:\t["
        << basepos
        << ":"
        << basepos+del_length-1
        << "]\n"
        << "Base pair length:\t"
        << del_length
        << endl
        << "Prev 10 base pairs:\t";
  printprevbases(finder1);
  cout  << endl
        << "region of interest:\t";
  while (printer != finder2 && printer != nullptr){
    cout  << printer->nuc;
    printer = printer->forward;
  }
  if (basepos < length && finder2 != nullptr) {
    finder2 = finder2->backward;
  }
  cout  << endl
        << "Next 10 base pairs:\t";
  if (basepos < length && finder2 != nullptr) {
    printnextbases(finder2);
  }
  else if (basepos == length){
    printnextbases(finder1);
  }
  else {
    printnextbases(b_end);
  }
  cout  << endl;
  pdeleter = finder1;
  if (basepos > 1 && basepos < length && finder2 != nullptr) {
    finder1 = finder1->backward;
    finder2 = finder2->forward;
    finder1->forward = finder2;
    finder2->backward = finder1;
    while (pdeleter != finder2) {
      base* next = pdeleter->forward;
      delete pdeleter;
      pdeleter = next;
    }
  }
  else if (basepos == 1) {
    finder1 = finder2->forward;
    b_start = finder1;
    finder1->backward = nullptr;
    finder2 = finder1;
    while (pdeleter != finder2) {
      base* next = pdeleter->forward;
      delete pdeleter;
      pdeleter = next;
    }
  }
  else if (basepos == length) {
    b_end = finder1->backward;
    b_end->forward = nullptr;
    delete pdeleter;
  }
  else if (finder2 == nullptr && basepos > 1 && basepos < length) {
    finder1 = finder1->backward;
    b_end = finder1;
    b_end->forward = nullptr;
    while (pdeleter != nullptr) {
      base* next = pdeleter->forward;
      delete pdeleter;
      pdeleter = next;
    }
  }
  else {
    cout  << "\nBASEPOS ERROR IN DELETION FUNCTION\n"
          << "CHECK INPUTVALIDATION\n";
  }
  cout << "\nThe above region of interest has been deleted\n";
  length = length-del_length;
}

// deletes a specified number of bases at a specified position and adds new bases
void DNAseq::replace(string baseinput, int basepos, int del_length, int choicemaker){
  deletion(basepos, del_length);
  add(baseinput, basepos, choicemaker);
  cout << "\nThe nucleotides have been replaced\n";
}

// save edited file from database to a fna file
void DNAseq::save(){
  string savefilename;
  base* seeker = b_start;
  bool inputvalid = true;
  string adder = ".fna";

  cout  <<"\nEnter a filename for the current DNA sequence to be saved in as a fna file:\n"
        << "> ";
  cin   >> savefilename;
  while (inputvalid) {
    if (savefilename.find(".") != -1) {
      cout  << "\nEnter filename again without formatting:\n"
            << "> ";
      cin   >> savefilename;
    }
    else {
      inputvalid = false;
      savefilename = savefilename + adder;
    }
  }
  ofstream tempstream (savefilename);
  tempstream  << ">File "
              << savefilename
              << endl;
  while (seeker != nullptr) {
    tempstream << seeker->nuc;
    seeker = seeker->forward;
  }
  cout << "\nThe DNA file has been saved\n";
}

// displays menu for sequence for editing
void DNAseq::submenu(){
  string userinput;
  string baseinput;
  string baseinput2;
  string baseinput3;
  string temp;
  int input = 0;
  int basepos;
  int basepos2;
  int basepos2a;

    cout  << "\nCURRENT DNA FILE: "
          << file_name
          << "\nSelect from one of the following options\n";
    while (input != 9){
      cout  << "1.\tFind DNA sequence by input\n"
            << "2.\tFind DNA sequence by file\n"
            << "3.\tAdd DNA sequence by input\n"
            << "4.\tAdd DNA sequence by file\n"
            << "5.\tDelete DNA sequence by input\n"
            << "6.\tReplace DNA sequence by input\n"
            << "7.\tReplace DNA sequence by file\n"
            << "8.\tSave edited DNA sequence\n"
            << "9.\tExit submenu\n"
            << "> ";
      cin   >> userinput;
      try {
        input = stoi(userinput);
        switch (input){
          case 1:{
            cout  << "\nEnter nucleotide sequence:\n"
                  << "> ";
            cin.ignore();
            getline(cin, baseinput);
            while (!validify(baseinput)){
              cout  << "\nPlease input a valid nucleotide sequence\n"
                    << "Valid bases are ACGTNUKSYMWRBDHV-\n"
                    << "> ";
              getline(cin, baseinput);
            }
            search(baseinput, 1); //WARNING CHECK INT VALUE
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 2:{
            cout  << "\nEnter file name nucleotide sequence:\n"
                  << "> ";
            cin   >> baseinput;
            temp = checker(baseinput);
            search(temp, 2);
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 3:{
            cout  << "\nEnter nucleotide sequence:\n"
                  << "> ";
            cin.ignore();
            getline(cin, baseinput);
            while (!validify(baseinput)){
              cout  << "\nPlease input a valid nucleotide sequence\n"
                    << "Valid bases are ACGTNUKSYMWRBDHV-\n"
                    << "> ";
              getline(cin, baseinput);
            }
            cout  << "\nEnter a base pair position:\n"
                  << "> ";
            getline(cin, baseinput2);
            basepos = stoi(baseinput2);
            while (validify2(basepos)) {
              cout  << "\nPlease input a valid base pair position\n"
                    << "Must be higher than 0 and lower than "
                    << length+2
                    << endl
                    << "> ";
              getline(cin, baseinput2);
              basepos = stoi(baseinput2);
            }
            add(baseinput, basepos, 1);
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 4:{
            cout  << "\nEnter file name nucleotide sequence:\n"
                  << "> ";
            cin   >> baseinput;
            temp = checker(baseinput);
            cout  << "\nEnter a base pair position:\n"
                  << "> ";
            getline(cin, baseinput2);
            basepos = stoi(baseinput2);
            while (validify2(basepos)) {
              cout  << "\nPlease input a valid base pair position\n"
                    << "Must be higher than 0 and lower than "
                    << length+2
                    << endl
                    << "> ";
              getline(cin, baseinput2);
              basepos = stoi(baseinput2);
            }
            add(temp, basepos, 2);
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 5:{
            cout  << "\nEnter a base pair position:\n"
                  << "> ";
            cin.ignore();
            getline(cin, baseinput);
            basepos = stoi(baseinput);
            while (validify(basepos)) {
              cout  << "\nPlease input a valid base pair position\n"
                    << "Must be higher than 0 and lower than "
                    << length+1
                    << endl
                    << "> ";
              getline(cin, baseinput);
              basepos = stoi(baseinput);
            }
            cout  << "\nEnter a base pair length:\n"
                  << "> ";
            getline(cin, baseinput2);
            basepos2 = stoi(baseinput2);
            basepos2a = validify3(basepos2, basepos);
            deletion(basepos, basepos2a);
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 6:{
            cout  << "\nEnter nucleotide sequence:\n"
                  << "> ";
            cin.ignore();
            getline(cin, baseinput);
            while (!validify(baseinput)){
              cout  << "\nPlease input a valid nucleotide sequence\n"
                    << "Valid bases are ACGTNUKSYMWRBDHV-\n"
                    << "> ";
              getline(cin, baseinput);
            }
            cout  << "\nEnter a base pair position:\n"
                  << "> ";
            getline(cin, baseinput2);
            basepos = stoi(baseinput2);
            while (validify(basepos)) {
              cout  << "\nPlease input a valid base pair position\n"
                    << "Must be higher than 0 and lower than "
                    << length+1
                    << endl
                    << "> ";
              getline(cin, baseinput2);
              basepos = stoi(baseinput2);
            }
            cout  << "\nEnter a base pair length:\n"
                  << "> ";
            getline(cin, baseinput3);
            basepos2 = stoi(baseinput3);
            basepos2a = validify3(basepos2, basepos);
            replace(baseinput, basepos, basepos2a, 1);
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 7:{
            cout  << "\nEnter file name nucleotide sequence:\n"
                  << "> ";
            cin   >> baseinput;
            temp = checker(baseinput);
            cout  << "\nEnter a base pair position:\n"
                  << "> ";
            getline(cin, baseinput2);
            basepos = stoi(baseinput2);
            while (validify(basepos)) {
              cout  << "\nPlease input a valid base pair position\n"
                    << "Must be higher than 0 and lower than "
                    << length+1
                    << endl
                    << "> ";
              getline(cin, baseinput2);
              basepos = stoi(baseinput2);
            }
            cout  << "\nEnter a base pair length:\n"
                  << "> ";
            getline(cin, baseinput3);
            basepos2 = stoi(baseinput3);
            basepos2a = validify3(basepos2, basepos);
            replace(temp, basepos, basepos2a, 2);
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 8:{
            save();
            system("PAUSE");
            cout  << endl;
          }
            break;
          case 9:{
            cout  << "Quitting...\n";
          }
            break;
          default:
            cout << "\nPlease use a valid input(3)\n";
            break;
        }
      }
        catch (exception &error) {
          cout  << "\nPlease use a valid input(4)\n"
                << error.what()
                << endl;
        }
    }
}

class DNADatabase{
  private:
    DNAseq* head;
    vector<string> filesowned;
    int filecount = 0;
  public:
    DNADatabase();
    int getfilecount();
    void option1();
    void option2();
    void option3(string specfile);
    void print();
    void load(string loadingDNAsequence);
};

// initialise pointer variable
DNADatabase::DNADatabase(){
  head = nullptr;
}

// get number of fies in database
int DNADatabase::getfilecount(){
  return filecount;
}

// dna_db points to a new sequence that points to the old sequence/null
void DNADatabase::load(string filename){
  head = new DNAseq(filename, head);
  filecount++;
}

// prints the name of each file in database
void DNADatabase::print(){
  for (int i=1; i<=filecount; i++){
    cout  << endl
          << "File "
          << i
          << ": "
          << filesowned[filecount-i];
  }
  cout << endl;
}

// search and validates for fna files then loads into database
// WARNING works with any number of files but only stores to max approx 82MB then unable to process when using more fna files, seems to be a problem with storage as chr1.fna (212MB) doesn't load
// get "bad alloc" error when dealing with large files
// WHY???????????
void DNADatabase::option1(){
  string user_in = "";
  string temp;
  string file;
  bool inputvalid = true;

  cin.ignore();
  while (inputvalid){
    getline(cin,user_in);
    // time.tic()                       //TIME
    istringstream iss(user_in);
    while(getline(iss,temp,' ')){
      istringstream iss2(temp);
      while(getline(iss2,file,',')){
        ifstream check(file, ios::in);
        if (check && (file.find(".fna")!=-1)) {
          load(file);
            inputvalid = false;
            filesowned.push_back(file);
            cout  << "\nFile "
                  << file
                  << " successfully loaded!";
        }
        else if ((!check) || (file.find(".fna")==-1)) {
          inputvalid = true;
          cout  << "\nUnable to open file:"
                << file
                << "\nPlease specify a valid file name with the correct format\n"
                << "> ";
        }
        else {
          cout << "OPTION 1 INPUT ERROR";
        }
      }
    }
  }
  // time.toc()                 //TIME
  cout  << endl
        << "\nNumber of files in database: "
        << filecount;
  /*cout  << endl
          << time
          << endl; */

}

// searches for specified sequence in database and accesses it for editing
void DNADatabase::option2(){
  string DNAinput = "0";
  int input = 0;
  bool inputvalid = true;

  cout << "\nSelect a DNA to process:\n";
  for (int i=1; i<=filecount; i++){
    cout  << i
          << ". "
          << filesowned[filecount-i]
          << endl;
  }
  cout  << "> ";
  cin   >> DNAinput;
  input = stoi(DNAinput);
  while (inputvalid){
    if (input<=filecount){
      inputvalid = false;
      if (input == 1){
        head->submenu();
      }
      else if (input > 1) {
        DNAseq* chooser = head->getnextseq();
        for (int i=1; i<(input-1); i++){
          chooser = chooser->getnextseq();
        }
        chooser->submenu();
      }
      else {
        cout << "\nERROR IN DNA SEQUENCE VALIDATER\n";
      }
    }
    else {
      cout  << "Please input a valid file number\n"
            << "> ";
      cin   >> DNAinput;
      input = stoi(DNAinput);
    }
  }
}

// searches entire database for matching sequence with specified file
void DNADatabase::option3(string specfile){
  string temp;
  ifstream seqanalyse;
  char c;
  char* tstart = nullptr;
  char* tpoint = nullptr;
  DNAseq* seeker;

  temp = checker(specfile);
  cout  << endl;
  seeker = head;
  for (int i=1; i<=filecount; i++){
    seeker->search(temp, 3); // WARNING CHECK INT VALUE
    seeker = seeker->getnextseq();
  }
}

int main(){
  DNADatabase dna_db;
  string userinput = "0";
  int input = 0;
  int counter = 0;
  string filename;

  /*          TIME
  TicToc time;
    Doesn't run????????
  */

  while (input != 4){
    system("CLS");
    cout  << endl
          << "Welcome to the DNA Editing program\n"
          << endl
          << "Select and option:\n"
          << "(1)\tLoad DNA(s).\n"
          << "(2)\tProcess a DNA.\n"
          << "(3)\tAnalyse the DNA database.\n"
          << "(4)\tQuit.\n"
          << "> ";
    cin   >> userinput;
    try {
      input = stoi(userinput);
      switch (input){
        case 1:{
          cout  << "\nEnter the DNA file names:\n"
                << "For multiple files, separate them by a comma. Only .fna are recognised.\n"
                << "> ";
          dna_db.option1();
          cout << endl;
          dna_db.print();
          cout  << endl;
          system("PAUSE");
        }
          break;
        case 2:{
          int loader = dna_db.getfilecount();
          if (loader == 0) {
            cout << "\nYou have not loaded in any files!\n";
          }
          else if (loader >= 1) {
            dna_db.option2();
          }
          else {
            cout << "ERROR IN CONTROL STATEMENT FOR DATABASE CHECK";
          }
          system("PAUSE");
          }
          break;
        case 3:{
            cout  << "Specify filename to be searched for inside the database:\n"
                  << "> ";
            cin   >> filename;
            // time.tic();            //TIME
            dna_db.option3(filename);
            // time.toc();            //TIME
            /*cout  << endl
                  << time
                  << endl;
            */
            system("PAUSE");
            }
          break;
        case 4:{
            cout  << "Quitting...\n";
            }
          break;
        default:{
        cout << "\nPlease use a valid input(1)\n";
        system("PAUSE");
      }
          break;
      }
    }
      catch (exception &error) {
        cout  << "\nPlease use a valid input(2)\n"
              << error.what();
        system("PAUSE");
      }
  }
  return 0;
}
