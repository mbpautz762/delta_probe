#include <iostream>
#include <fstream>
using namespace std;

int N, N_delta, K = 100;
string *sequence;
bool *is_delta;

struct probe_info {
  int seq_index, pos;
  int fp_count, fn_count;
  double error_rate;
};

void read_input(void)
{
  ifstream input("covid_small.txt");
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

// Does position p in sequence[i] match position q in sequence[j] (1 mismatch char ok)?
bool is_match(int i, int p, int j, int q)
{
  int mismatched = 0;
  for (int k=0; k<K; k++)
    if (sequence[i][p+k] != sequence[j][q+k]) {
      mismatched ++;
      if (mismatched > 1) return false;
    }
  return true;
}

probe_info eval_probe(int i, int p)
{
  probe_info info;
  info.seq_index = i;
  info.pos = p;
  info.fp_count = info.fn_count = 0;

  for (int j=0; j<N; j++) {
    if (is_delta[j]) info.fn_count ++;
    for (int q=0; q<sequence[j].length()-K; q++) {
      if (is_match(i,p,j,q)) {
	if (is_delta[j]) info.fn_count --;
	if (!is_delta[j]) info.fp_count ++;
	break;
      }
    }
  }

  double FPR = (double)info.fp_count / (N-N_delta);
  double FNR = (double)info.fn_count / N_delta;
  info.error_rate = 2.0 * FPR + 1.0 * FNR;
  return info;
}

int main(void)
{
  read_input();

  // Loop over all possible probes, remember the best one...
  probe_info best;
  best.error_rate = 99999;
  for (int i=0; i<N; i++) {
    cerr << ".";
    if (is_delta[i])
      for (int p=0; p<sequence[i].length()-K; p++) {
	// Evaluate probe at position p in sequence[i]
	probe_info info = eval_probe(i,p);
	if (info.error_rate < best.error_rate) best = info;
      }
  }
  cerr << "\n";
  
  // Print out info about best
  cout << "Best probe: " << sequence[best.seq_index].substr(best.pos, K) << "\n";
  cout << "False positives: " << best.fp_count << "\n";
  cout << "False negatives: " << best.fn_count << "\n";
  cout << "Error_rate: " << best.error_rate << "\n";
  cout << "Worst-case runtime: Theta(N^2 M^2 K)\n";   // Runtime in terms of N, M, K
  cout << "Anticipated runtime: Theta(N^2 M^2)\n";   // Runtime in terms of N, M, K
  
  delete [] sequence;
  delete [] is_delta;
}
