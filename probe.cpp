#include <iostream>
#include <fstream>
using namespace std;

// M is number of buckets, and number to modulo by during hashing
// X is the prime constant used in hashing functions.
// xTerm is set to initial value needed to pre-computer X^K-1 in calcPower
const int K = 50;
const int M = 1000000007;
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
    // simplified constructor for only holding the sequence info
    pNode(int seq, pNode *n) : sequence(seq), next(n) {;}

};

struct Node {
    int key;
    pNode *locations;
    Node *next;

    Node(int k, int seq, int pos, Node *n) : key(k), next(n) { locations = new pNode(seq, pos, nullptr); }
};

struct Probe {
    int lH;
    int rH;
    int firstSeq;
    int firstPos;
    int falsePos = 0;
    int falseNeg = N_delta;
    double FPR;
    double FNR;
    double errorRate;
    bool* seqMatches;

    Probe(int n) {
        seqMatches = new bool[n];
    }
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

void updateProbeInfo(Probe *ptr) {
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
// checks if a key is in the hash table.  if it is, that key's
// location info is added into it's pNode.  if not found, 
// a new node is added
bool updateTable(Node **table, int h, int i, int j) {
    if (table[h] == nullptr) {
        table[h] = new Node(h, i, j, table[h]);
        return true;
    }

    else {
        Node *tmp = table[h];
        while (tmp) {
            if (tmp->key == h) {
                tmp->locations = new pNode(i, j, tmp->locations);
                return true;
            }
            else tmp = tmp->next;
        }
        // if matching key not already in table, it is a collision.
        // add a new link in the node
        table[h] = new Node(h, i, j, table[h]);
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

int hashChecker(string *seq, int hash, int pos) {
    int h = 0;
    for (int i = 0; i < K; i++) {
        h = ((((long long)h * X) % M) + seq->at(i + pos)) % M;
    }

    return h;
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
                // cout << "rolling: " << h;
                // cout << "  checked: " << hashChecker(&sequence[i], h, j) << endl;
                
                updateTable(table, h, i, j);

            }            
                
    }

    int lHash, rHash;
    Probe *currentProbe = new Probe(N);
    Probe *bestProbe = new Probe(N);
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
                rHash = hasher(&sequence[i], rHash, (j + K));
                // cout << "lhash: " << lHash << " checked: " << hashChecker(&sequence[i], lHash, j) << endl;
                // cout << "rhash: " << rHash << " checked: " << hashChecker(&sequence[i], rHash, (j + K)) << endl;

                currentProbe->firstSeq = i;
                currentProbe->firstPos = j;
                currentProbe->lH = lHash;
                currentProbe->rH = rHash;
                // first, check all left hash locations for a matching right hash
                if (table[lHash]) {
                    Node *match = table[lHash];
                    // advance to the node of the key you're checking (in case of collisions)
                    while (match->key != lHash) {
                        match = match->next;
                    }

                    pNode *locs = match->locations;
                    while (locs) {
                        // NOTE: I reversed the order of the arguments in isMatch and for some reason it works better than before.  WHY
                        if (isMatch(locs->sequence, locs->position + K, i, j)) {
                            // isMatch(i, j + K, locs->sequence, locs->position + K)) {
                            // if match found in delta location, decrement falsePos, 
                            // otherwise, increment falseNeg.  Then add location to 
                            // bool array to track when checking the right hashes
                            if (is_delta[locs->sequence]) (currentProbe->falseNeg)--;
                            else (currentProbe->falsePos)++;
                            currentProbe->seqMatches[locs->sequence] = true;
                        }

                        locs = locs->next;
                    }
                }
                // now, check all the right hashes.  only update statistics if a match was
                // not already found on the sequence
                if (table[rHash]) {
                    Node *match = table[rHash];
                    // advance to the node of the key you're checking (in case of collisions)
                    while (match->key != rHash) {
                        match = match->next;
                    }
                    pNode *locs = match->locations;
                    while (locs) {
                        // skip locations on sequences that have already matched
                        // take advantage of short circuiting here
                        if (currentProbe->seqMatches[locs->sequence] == false && isMatch(locs->sequence, locs->position - K, i, j - K)) {
                        // isMatch(i, j - K, locs->sequence, locs->position - K)) {
                            if (is_delta[locs->sequence]) (currentProbe->falseNeg)--;
                            else (currentProbe->falsePos)++;
                            currentProbe->seqMatches[locs->sequence] = true;
                        }

                        locs = locs->next;                   

                    }
                }
                updateProbeInfo(currentProbe);
                if (currentProbe->errorRate < bestProbe->errorRate) bestProbe = currentProbe;
                
                   




            
            }            
        }

    }

    cout << "Best Node: " << endl;
    cout << "   probe:      " << sequence[bestProbe->firstSeq].substr(bestProbe->firstPos, K * 2) << endl;
    cout << "   left hash:  " << bestProbe->lH << endl;
    cout << "   left hash:  " << bestProbe->rH << endl;
    cout << "   false pos:  " << bestProbe->falsePos << endl;
    cout << "   false neg:  " << bestProbe->falseNeg << endl;
    cout << "   error rate: " << bestProbe->errorRate << endl;


    deleteTable(table);

    return 0;

}