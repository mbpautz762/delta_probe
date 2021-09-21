#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// M is number of buckets, and number to modulo by during hashing
// X is the prime constant used in hashing functions.
// xTerm is set to initial value needed to pre-computer X^K-1 in calcPower
const int K = 50;
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
    bool evaluated = false;
    int falsePos = 0;
    int falseNeg = N_delta;
    double FPR;
    double FNR;
    double errorRate;
    pNode *locations;
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

int calcPower(int length) {
    for (int i = 0; i < length - 1; i++) {
        xTerm = ((long long)xTerm * X) % M;
    }

    return xTerm;
}

bool updateTable(Node **table, int &h, int &i, int &j) {
    // if node doesn't exist, create a new one.
    if (table[h] == nullptr) {
        table[h] = new Node(h, i, j, table[h]);
    }
    // if hashed index exists, instead of making a new node, 
    // just add a new link in the pNode list to add the additional position information
    else if (table[h]) {
        table[h]->locations = new pNode(i, j, table[h]->locations);
    }

    return true;

}

int hasher(string *seq, int hash, int pos) {
    // if passing in a 0 hash, calculate the entire initial hash
    // also, calculate x^N - 1
    if (hash == 0) {
        for (int i = 0; i < K; i++) {
            hash = ((((long long)hash * X) % M) + seq->at(i)) % M;
            //DEBUG: checking to make sure hash doesn't get screwy
            if (hash <= 0){
                cout << "Error! hash is not calculating properly." << endl;
                exit (1);
            }
        }
    }
    
    // if hash already exists, calculate the rolling hash instead
    // deletes and adds relative to pos - 1 to account for the initialization above
    else {
        hash = (((hash + M - (((long long)seq->at(pos - 1) * xTerm) % M) % M) * X) % M
            + seq->at((pos - 1) + K)) % M;
    }

    return hash;
}
// Does position p in sequence[i] match position q in sequence[j] (1 mismatch char ok)?
bool isMatch(int sourceSeq, int sourcePos, int compareSeq, int comparePos) {
      int mismatched = 0;
  for (int k=0; k<K; k++)
    if (sequence[sourceSeq][sourcePos+k] != sequence[compareSeq][comparePos+k]) {
      mismatched++;
      if (mismatched > 1) return false;
    }
  return true;
}


int main(void) {
    readInput();

    Node **table = new Node *[M] {nullptr};

    int i = 0, j = 0, h = 0;
    // pre-calc x^N - 1
    xTerm = calcPower(K);

    // loop through all sequences, delta and covid, hashing every position in every sequence
    for (i = 0; i < N; i++) {
        // hash each position of the sequence to table, storing 
        // sequence and position data in pNode.
            h = 0;
            for (j = 0; ((j + K) < sequence[i].length()); j++) {
                h = hasher(&sequence[i], h, j);
                
                updateTable(table, h, i, j);

                // table[h]->falseNeg--;
                // updateProbeInfo(table[h]);
                // if (table[h]->errorRate < bestProbe->errorRate) bestProbe = table[h];
            }            
                
    }
    // now, test all covid sequences to check for false positives
    // for (i = 0; i < N; i++) {
    //     if (is_delta[i] == false) {
    //         h = 0;
    //         for (j = 0; (j + K) < sequence[i].length(); j++) {
    //             h = hasher(&sequence[i], h, j);
    //             // don't actually create nodes.  only check if
    //             // it hashes to an existing delta hash
    //             if (table[h]) {
    //                 table[h]->falsePos++;

    //                 updateProbeInfo(table[h]);
    //                 if (table[h]->errorRate < bestProbe->errorRate) bestProbe = table[h];

    //             }   
    //         }   
    //     }
    // }
    int lHash, rHash;
    Node *currentProbe = new Node(0, 0, 0, nullptr);
    Node *bestProbe = new Node(0, 0, 0, nullptr);
    bestProbe->errorRate = 1000;
    Node *probe;
    pNode *locations;   


    for (int i = 0; i < N; i++) {
        // check every position in every delta sequence
        if (is_delta[i]) {
            // start with a fresh hash for hashing function
            lHash = 0, rHash = 0;
            for (int j = 0; j + (2 * K) < sequence[i].length(); j++) {

                lHash = hasher(&sequence[i], lHash, j);
                rHash = hasher(&sequence[i], rHash, j);
                
                // ONLY evaluate if that position has not been checked elsewhere
                // first, look at lHash and all its matching right side locations
                if (table[lHash]->evaluated == false) {
                    probe = table[lHash];
                    locations = probe->locations;
                    // TODO: Fix it to where it stops looking for matches if a sequence has already been matched
                    while (locations) {
                        if (isMatch(i, j + K, locations->sequence, locations->position + K)) {
                            int matchingSequence = locations->sequence;
                            // if matching location belongs to a delta sequence, decrement falseNeg
                            // otherwise, increment false pos
                            if (is_delta[locations->sequence]) (probe->falseNeg)--;
                            else (probe->falsePos)++;

                            // if match found, fast forward to the next sequence
                            while (locations && locations->sequence == matchingSequence) locations = locations->next;
                            
                        }
                        else locations = locations->next;
                    }
                    // after looping through all the locations, update stats and compare with bestNode

                }
                    updateProbeInfo(probe);
                    if (probe->errorRate < bestProbe->errorRate) bestProbe = probe;                    
                    probe->evaluated = true;




            
            }            
        }

    }

    cout << "Best Node: " << endl;
    cout << "   key:        " << bestProbe->key << endl;
    cout << "   probe:      " << sequence[bestProbe->locations->sequence].substr(bestProbe->locations->position, K * 2) << endl;
    cout << "   false pos:  " << bestProbe->falsePos << endl;
    cout << "   false neg:  " << bestProbe->falseNeg << endl;
    cout << "   error rate: " << bestProbe->errorRate << endl;

    cout << "   all locations:  " << endl;
    pNode *tmp = bestProbe->locations;
    while (tmp) {
        cout << "sequence:  " << tmp->sequence << " position:   " << tmp->position << endl;
        tmp = tmp->next;
    }

    deleteTable(table);

    return 0;

}