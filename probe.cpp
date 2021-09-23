#include <iostream>
#include <fstream>
using namespace std;

// global variables
// P is prime number to modulo by during hashing
// X is the prime constant used in hashing functions
// xTerm is set to initial value needed to pre-compute X^K-1 in calcPower
// K is the probe length
const int K = 100;
const int P = 1000000007;
const int X = 5;
int xTerm = 1;
int hashSize = K / 2;

int N, N_delta;
string *sequence;
bool *is_delta;

void readInput(void)
{
  ifstream input("covid.txt");
//   ifstream input("covidsmall.txt");
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
// holds location information for Node
struct pNode {
    int sequence;
    int position;
    pNode *next;

    pNode(int seq, int pos, pNode *n) : sequence(seq), position(pos), next(n) {;}
    // simplified constructor for only holding the sequence info
    pNode(int seq, pNode *n) : sequence(seq), next(n) {;}

};
// holds a hashed sequence, and its associated location information
struct Node {
    int key;
    pNode *locations;
    Node *next;

    Node(int k, int seq, int pos, Node *n) : key(k), next(n) { locations = new pNode(seq, pos, nullptr); }
    Node(int k, pNode *p, Node *n) : key(k), next(n), locations(p) {;}
};
class Table {
    private:
        int N = 0;
        // initial hash table size
        int M = 50000;     

        Node **table;

    public:
        Table();
        ~Table();

        Node *find(int key);
        void insert(int key, int seq, int pos);
        void print(void);
        int hash(int key, int tableSize) { return key % tableSize; }
        Node **getTable() { return table; }
        Node *getTableEntry(int key) { return table[hash(key, M)]; }
        int getN() { return N; }
        void deleteTable();
};

Table::Table() { table = new Node *[M] {nullptr}; }
Table::~Table() { deleteTable(); }

// search for a key in a hashed bucket
// if found, returns a pointer to the containing node
Node *Table::find(int key) {
    int h = hash(key, M);
    if (table[h]) {
        Node *tmp = table[h];
        while (tmp) {
            if (key == tmp->key) return tmp;
            tmp = tmp->next;
        }
    }
    // if no match, return null
  return nullptr;
}

void Table::insert(int key, int seq, int pos) {
    int i, h;

    // check if the table needs to be upsized at time of insert    
    if (N > (2 * M) ) {
        Node **newTable = new Node *[M * 2] {nullptr};
        Node *tmp;
        Node *prev;

        for (int i = 0; i < M; i++) {
            tmp = table[i];
            while (tmp) {
                h = hash(tmp->key, (M * 2));
                newTable[h] = new Node(tmp->key, tmp->locations, newTable[h]);
                // while advancing to next node, free the memory
                // at current position in the old table.
                // since we're already here, why not?
                prev = tmp;
                tmp = tmp->next;
                delete prev;            }                  
        }
        // double the table size AFTER re-hashing everything
        M *= 2;

        delete [] table;
        table = newTable;
    }
    // now we insert the new node.
    // if key already exists in table, just add new location info
    N++;

    h = hash(key, M);
    if (table[h] == nullptr || (!find(key)) ) {
        table[h] = new Node(key, seq, pos, table[h]);
    }
    else {
        Node *tmp = find(key);
        tmp->locations = new pNode(seq, pos, tmp->locations);
    }
}
// more for debugging than anything else
void Table::print(void) {
    Node *tmp;
    pNode *p;
        for (int i=0; i<M; i++){
        if (table[i]) {
            tmp = table[i];
            p = table[i]->locations;
            while (tmp) {
                cout << "\nkey: " << tmp->key;
            while (p) {
                cout << " locs: " << p->sequence << " " << p->position;
                p = p->next;
            }
            tmp = tmp->next;
            }
        }
    }
}
// systematically frees all nodes, and all associated pNodes
void Table::deleteTable() {
    Node *current;
    Node *prev;
    pNode *loc;
    pNode *locPrev;
    for (int i = 0; i < M; i++) {
        current = table[i];
        while (current) {
            prev = current;
            if (current->locations) loc = current->locations; {
                while (loc) {
                    locPrev = loc;
                    loc = loc->next;
                    delete locPrev;
                }
            }

            current = current->next;
            delete prev;
        }            
    }
    delete[] table;
    return;
}
// holds information about the K-length probe
class Probe {
    public:
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

        Probe() { seqMatches = new bool[N] {false}; }
        ~Probe() { ; }

        Probe(int lhash, int rhash, int i, int j) : lH(lhash), rH(rhash), firstSeq(i), firstPos(j) {seqMatches = new bool[N] {false};}
        // int getlHash() { return lH; }
        // int getrHash() { return rH; }
        // int getSeq() { return firstSeq; }
        // int getPos() { return firstPos; }
        // int getFalsePositives() { return falsePos; }
        // int getFalseNegatives() { return falseNeg; }
        // double getFPR() { return FPR; }
        // double getFNR() { return FNR; }
        // double getErrorRate() { return errorRate; }
        // bool getSeqMatches() { return seqMatches; }
        // bool is_SeqMatched(int i) { return seqMatches[i]; }

        // // setters
        // void setSeqMatch(int i) { seqMatches[i] = true;}
        // void setlHash(int h) { lH = h; }
        // void setrHash(int h) { rH = h; }
        // void setSeq(int seq) { firstSeq = seq; }
        // void setPos(int pos) { firstPos = pos; }
        // void decrementFN() { falseNeg--; }
        // void incrementFP() { falsePos++ ;}

        bool evaluate(Node *match, int seq, int pos);
        void updateInfo();
        void reset(int lHash, int rHash, int i, int j); 
        void cleanup();    
};

// evaluates a probe's effectiveness
// tracks match locations to prevent double matching
bool Probe::evaluate(Node *match, int seq, int pos) {
    bool isMatch(int sourceSeq, int sourcePos, int compareSeq, int comparePos);
    pNode *locs = match->locations;
    while (locs) {
        if (seqMatches[locs->sequence] == false 
            && isMatch(seq, pos, locs->sequence, pos)) {
            // if match found in delta location, decrement falsePos, 
            // otherwise, increment falseNeg.  Then add location to 
            // bool array to track when checking the right hashes
            if (is_delta[locs->sequence]) (falseNeg)--;
            if (!is_delta[locs->sequence]) (falsePos)++;

            seqMatches[locs->sequence] = true;
        }
        locs = locs->next;
    }        
    return true;
}

void Probe::updateInfo() {
    FPR = ((double)(falsePos) / (N - N_delta));
    FNR = ((double)(falseNeg) / (N_delta));
    errorRate = (2.0 * FPR) + (1.0 * FNR);

    return;
}

int calcPower(int length) {
    for (int i = 0; i < length - 1; i++) {
        xTerm = ((long long)xTerm * X) % P;
    }

    return xTerm;
}

int hasher(string *seq, int hash, int pos) {
    // if passing in a 0 hash, calculate the entire initial hash
    // also, calculate x^N - 1
    if (hash == 0) {
        for (int i = 0; i < hashSize; i++) {
            hash = ((((long long)hash * X) % P) + seq->at(pos + i)) % P;
            //DEBUG: checking to make sure hash doesn't get screwy
        }
    }
    
    // if hash already exists, calculate the rolling hash instead
    // deletes and adds relative to pos - 1 to account for the initialization above
    else {
        hash = (((hash + P - (((long long)seq->at(pos - 1) * xTerm) % P) % P) * X) % P
            + seq->at((pos - 1) + (hashSize))) % P;        
    }

    return hash;
}
// re-initializes probe with fresh info
void Probe::reset(int lHash, int rHash, int i, int j) {
    lH = lHash;
    rH = rHash;
    firstSeq = i;
    firstPos = j;
    falsePos = 0;
    falseNeg = N_delta;
    delete [] seqMatches;
    seqMatches = new bool[N] {false};
}

// clears bool array and frees memory.  used instead of destructor to preven
// a double-free issue
void Probe::cleanup(){
    delete [] seqMatches;
}


// The expensive check: checks a matching location char by char for a match,
// allowing up to one mismatching character.
bool isMatch(int sourceSeq, int sourcePos, int compareSeq, int comparePos) {
    int mismatched = 0;
    for (int k=0; k < hashSize; k++) {
        if (sequence[sourceSeq][sourcePos+k] != sequence[compareSeq][comparePos+k]) {
            mismatched++;
            if (mismatched > 1) return false;
        }
    }
    return true;
}


int main(void) {

    readInput();

    Table table;
    int i = 0, j = 0, h = 0;

    // pre-calc x^N - 1 for hashing function
    xTerm = calcPower(hashSize);

    // loop through all sequences, delta and covid, hashing each position in every sequence.  
    // sequence and position information stored in pNode
    for (i = 0; i < N; i++) {
        h = 0;
        for (j = 0; ((j + hashSize) < sequence[i].length()); j++) {
            h = hasher(&sequence[i], h, j);
            table.insert(h, i, j);
        }            
                
    }
    int lHash, rHash;
    Probe bestProbe;
    Probe currentProbe;
    bestProbe.errorRate = 9999;

    for (int i = 0; i < N; i++) {
        // check every position in every delta sequence
        if (is_delta[i]) {
            // start with a fresh hash - required for hashing function
            lHash = 0, rHash = 0;
            for (int j = 0; j + K < sequence[i].length(); j++) {
                lHash = hasher(&sequence[i], lHash, j);
                rHash = hasher(&sequence[i], rHash, (j + hashSize));
                currentProbe.reset(lHash, rHash, i, j);

                // first, check all left hash locations for a matching right hash
                Node *match = table.getTableEntry(lHash);

                // advance to the node of the key you're checking (in case of collisions)
                while (match->key != lHash) {
                    match = match->next;
                }
                currentProbe.evaluate(match, i, j + hashSize);

                // now, check all the right hashes.  only update statistics if a match was
                // not already found on that particular sequence
                match = table.getTableEntry(rHash);
                while (match->key != rHash) {
                    match = match->next;
                }
                currentProbe.evaluate(match, i, j);

                // all evaluations are done - update probe statistics and compare
                currentProbe.updateInfo();

                if (currentProbe.errorRate < bestProbe.errorRate) bestProbe = currentProbe;
            }            
        }
    }
    cout << "Best Probe: " << endl;
    cout << "   probe:      " << sequence[bestProbe.firstSeq].substr(bestProbe.firstPos, K) << endl;
    cout << "   left hash:  " << bestProbe.lH << endl;
    cout << "   left hash:  " << bestProbe.rH << endl;
    cout << "   false pos:  " << bestProbe.falsePos << endl;
    cout << "   false neg:  " << bestProbe.falseNeg << endl;
    cout << "   error rate: " << bestProbe.errorRate << endl;
    cout << "   sequence:   " << bestProbe.firstSeq << endl;
    cout << "   position:   " << bestProbe.firstPos << endl;

    delete[] sequence;
    delete[] is_delta;

    currentProbe.cleanup();

    return 0;

}