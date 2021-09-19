#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// 10^9 + 9 originally chosen for M based on research on best modulo to prevent collisions.
// M is number of buckets, and number to modulo by during hashing
// X is the prime constant used in hashing functions.
// xTerm is initial value used to pre-computer X^K-1 
const int K = 100;
const int M = 100000007;
const int X = 5;
int xTerm = 1;

int N, N_delta;
string *sequence;
bool *is_delta;

struct pNode {
    int sequence;
    int position;
    pNode *next;

    pNode(int seq, int pos, pNode *n) : sequence(seq), position(pos), next(n) {;}

};

struct Node {
    int key;
    int falsePos = 0;
    int falseNeg = N_delta;
    double FPR;
    double FNR;
    double errorRate;
    pNode * locations;
    Node *next;

    Node(int k, int seq, int pos, Node *n) : key(k), next(n) { locations = new pNode(seq, pos, nullptr); }
};

void deleteTable(Node **table) {
    Node *current;
    Node *prev;
    pNode *pPtr;
    // TODO: figure out how to delete pNode memory properly
    // and add to this funciton
    for (int i = 0; i < M; i++) {
        current = table[i];
        if (current) {
            prev = current;
            pPtr = current->locations;
            current = current->next;
            delete prev;
            delete pPtr;
        }
    }
    delete[] table;

    return;
}

void updateProbeInfo(Node *ptr) {
    ptr->FPR = ((double)(ptr->falsePos) / (N - N_delta));
    ptr->FNR = ((double)(ptr->falseNeg) / (N_delta));
    ptr->errorRate = 2.0 * (ptr->FPR) + 1.0 * (ptr->FNR);
}

void readInput(void)
{
  ifstream input("covid.txt");
  string label, seq;
  while (input >> label >> seq) N++;
  input.clear();
  input.seekg(0);
  sequence = new string[N];
  is_delta = new bool[N];
  for (int i=0; i<N; i++) {
    input >> label >> sequence[i];
    is_delta[i] = label == "delta_variant";
    if (is_delta[i]) N_delta++;
  }
}

int strHash(string key) {
    unsigned long h = 0;

    // Horner's Rule:
    for (int i = 0; i < key.length(); i++) {
        h = ((((long long)h * X) % M) + key.at(i)) % M;        
    }
    return h;
}

int calcPower(int length) {
    for (int i = 0; i < length - 1; i++) {
        xTerm = ((long long)xTerm * X) % M;
    }

    return xTerm;
}

int rollingHash(string key, int hash, int pos, int length) {
    if (xTerm == 1) cout << "Warning!  X^N - 1 is not pre-computed" << endl;
    
    int h = hash;
        h = (((h + M - (((long long)key.at(pos) * xTerm) % M) % M) * X) % M
            + key.at(pos + length)) % M;

    return h;

}

int main(void) {
    readInput();

    Node **table = new Node *[M] {nullptr};
    Node *bestNode = new Node(0, 0, 0, nullptr);   

    int i = 0, j = 0, h = 0;

    // loop through all delta sequences
    for (i = 0; i < N; i++) {
        // hash a substring of the sequence to table, storing 
        // sequence and position data in pNode.
        if (is_delta[i]) {

            for (j = 0; (j + K) <= sequence[i].length(); j++) {
                h = strHash(sequence[i].substr(j, K));
                // if node doesn't exist, create a new one.
                if (table[h] == nullptr) {
                    table[h] = new Node(h, i, j, table[h]);
                }
                // if hashed index exists, instead of making a new node, 
                // just add a link to the pNode to add the additional position information
                else if (table[h]) {
                    table[h]->locations = new pNode(i, j, table[h]->locations);
                }
                // either way, this hashed location is a 'match', so its associated
                // false negative value should be decremented.
                table[h]->falseNeg--;
                updateProbeInfo(table[h]);
                if (table[h]->errorRate < bestNode->errorRate) bestNode = table[h];
            }            
        }        
    }
    // now, test all covid sequences to check for false positives
    for (i = 0; i < N; i++) {
        if (is_delta[i] == false) {
            for (j = 0; (j + K) <= sequence[i].length(); j++) {
                h = strHash(sequence[i].substr(j, K));
                // don't actually create nodes.  only check if
                // it hashes to an existing delta hash
                if (table[h]) {
                    table[h]->falsePos++;

                    updateProbeInfo(table[h]);
                    if (table[h]->errorRate < bestNode->errorRate) bestNode = table[h];

                }   
            }   
        }
    }

    cout << "Best Node: " << endl;
    cout << "   key:        " << bestNode->key << endl;
    cout << "   probe:      " << sequence[bestNode->locations->sequence].substr(bestNode->locations->position, K) << endl;
    cout << "   false pos:  " << bestNode->falsePos << endl;
    cout << "   false neg:  " << bestNode->falseNeg << endl;
    cout << "   error rate: " << bestNode->errorRate << endl;

    cout << "   all locations:  " << endl;
    pNode *tmp = bestNode->locations;
    while (tmp) {
        cout << "sequence:  " << tmp->sequence << " position:   " << tmp->position << endl;
        tmp = tmp->next;
    }

    deleteTable(table);

    return 0;

}