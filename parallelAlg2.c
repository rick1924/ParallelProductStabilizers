/*Function that performs the H-CX-CZ-S-H circuit extraction (Algorithm 2 of the paper) and reduction to basis state

Constraints:
 - We have to make sure that N>=M
 - num_qubits % N = 0

 Inputs:
 - A stabilizer tableau in its canonical reduced form
 - A vector of phases corresponding to the rows of the tableau

 Outputs:
 - A csv file containing the matrix corresponding to a basis state and its vector of phases
*/ 

#include "stdio.h"
#include "stdlib.h"
#include "bsp.h"
#include "math.h"
#include "string.h"
#include "stdbool.h"

// Modulo function that works for negative numbers
#define MOD(a, b) (((a % b) + b) % b)

const int P = 4;  // Number of processors
const int N = 2;  // Processor rows
const int M = 2;  // Processor columns
const int num_qubits = 8; // Sometimes referred to as n
int our_nan = -1 * ((2 * num_qubits) + 1);
bool has_one[2 * 8];
int mat[8][2 * 8];  // Size of the tableau
int vec[8];  // Size of the phase vector
bool found_a_one = false;

// Define a structure used to represent a single-qubit or two-qubit gates
struct gate{
  char name[3];
  int i;
  int j;
};
// Counters for the number of gates we have of each kind
int h1counter = 0;
int cxcounter = 0;
int czcounter = 0;
int scounter = 0;
int h2counter = 0;

struct gate had1[num_qubits * num_qubits];
struct gate cnot[num_qubits * num_qubits];
struct gate cphase[num_qubits * num_qubits];
struct gate phase[num_qubits];
struct gate had2[num_qubits];

// Functions to append gates to a list
void h1append(int i){
  strcpy(had1[h1counter].name, "H");
  had1[h1counter].i = i;
  had1[h1counter].j = -1;
  h1counter++;
}

void h2append(int i){
  strcpy(had2[h2counter].name, "H");
  had2[h2counter].i = i;
  had2[h2counter].j = -1;
  h2counter++;
}

void cxappend(int i, int j){
  strcpy(cnot[cxcounter].name, "CX");
  cnot[cxcounter].i = i;
  cnot[cxcounter].j = j;
  cxcounter++;
}

void czappend(int i, int j){
  strcpy(cphase[czcounter].name, "CZ");
  cphase[czcounter].i = i;
  cphase[czcounter].j = j;
  czcounter++;
}

void sappend(int i){
  strcpy(phase[scounter].name, "S");
  phase[scounter].i = i;
  phase[scounter].j = -1;
  scounter++;
}

// Function that reads the tableau in its canonical form and loads it into directly into our 2-d array "mat"
void read_mat(){

  char buffer[100000];
  char *record, *line;
  int i = 0, j = 0;
  FILE *fstream = fopen("~/ParallelStabilizerInnerProduct/Outputs/alg1_output_mat.txt", "r");
  if (fstream == NULL){
    printf("\n file opening failed ");
  }
  while ((line = fgets(buffer, 100000, fstream)) != NULL){
    record = strtok(line, ",");
    while (record != NULL){
      mat[i][j++] = atoi(record);
      record = strtok(NULL, ",");
    }
    ++i;
  }

  fclose(fstream);
}

// Function that loads the vector of phases directly into "vec"
void read_vec(){

  char buffer[100000];
  char *record, *line;
  int i = 0;
  FILE *fstream = fopen("~/ParallelStabilizerInnerProduct/Outputs/alg1_output_vec.txt", "r");
  if (fstream == NULL){
    printf("\n file opening failed ");
  }
  while ((line = fgets(buffer, 100000, fstream)) != NULL){
    record = strtok(line, ",");
    while (record != NULL){
      vec[i++] = atoi(record);
      record = strtok(NULL, ",");
    }
  }

  fclose(fstream);
}

// Function that writes the tableau and resulting phases into a CSV file
void write_output(){
  FILE *alg2_output_mat;
  FILE *alg2_output_vec;

  alg2_output_mat = fopen("~/ParallelStabilizerInnerProduct/alg2_output_mat.txt", "w");

  if (alg2_output_mat != NULL){
    for (int i = 0; i < num_qubits; i++){
      for (int j = 0; j < 2 * num_qubits; j++){
        if (i == num_qubits - 1 && j == 2 * num_qubits - 1){
          fprintf(alg2_output_mat, "%d \n", mat[i][j]);
        }
        else{
          fprintf(alg2_output_mat, "%d,", mat[i][j]);
        }
      }
    }
  }
  else{
    printf("The file could not be opened");
  }

  fclose(alg2_output_mat);
  alg2_output_vec = fopen("~/ParallelStabilizerInnerProduct/Outputs/alg2_output_vec.txt", "w");

  if (alg2_output_vec != NULL){
    for (int l = 0; l < num_qubits; l++){
      if (l != num_qubits - 1){
        fprintf(alg2_output_vec, "%d,", vec[l]);
      }
      else{
        fprintf(alg2_output_vec, "%d", vec[l]);
      }
    }
  }

  fclose(alg2_output_vec);
}

// Parallel function
void parallelalg2()
{

  // Begin using "P" processors
  bsp_begin(P);

  // Name the processors using 2d naming 
  int pid = bsp_pid(); // P0 = 00, P1 = 10, P2 = 01, P3 = 11
  int row_name = pid % N;
  int col_name = (int)floor(pid / N);

  const int column_dim = 2 * num_qubits / M;
  const int row_dim = num_qubits / N;
  const int tot_values = column_dim * row_dim;

  //  Initialize "my_matrix"
  int *my_matrix[row_dim];
  for (int i = 0; i < row_dim; i++){
    my_matrix[i] = malloc(column_dim * sizeof(int));
  }

  // Distribute the data of the global "mat" to the local matrices "my_matrix" according to the 2d cyclic Cartesian distribution
  for (int i = 0; i < row_dim; i++){
    for (int j = 0; j < column_dim; j++){
      my_matrix[i][j] = mat[row_name + i * N][col_name + j * M];
    }
  }

  // Initialize "my_vec" and fill it with entries from "vec"
  int *my_vec = NULL;
  if (col_name == 0){
    my_vec = malloc(row_dim * sizeof(int));
    for (int i = 0; i < row_dim; i++){
      my_vec[i] = vec[row_name + i * N];
    }
  }

  // Initialize registers that enable data communication between the processors
  int *entries_received = malloc(column_dim * sizeof(int));
  bsp_push_reg(entries_received, column_dim * sizeof(int)); // Array for received row entries. Used in swapping
  int phase_received = -1;
  bsp_push_reg(&phase_received, sizeof(int));  // Array for received phase entries

  int *winner_row_array = malloc(P * sizeof(int));
  for (int i = 0; i < P; i++){
    winner_row_array[i] = -1;
  }
  bsp_push_reg(winner_row_array, P * sizeof(int)); // Array to determine which processor wins
  int winner_row = -1;
  bsp_push_reg(&winner_row, sizeof(int)); // Which row wins
  int *first_row_array = malloc(num_qubits * sizeof(int));
  for (int i = 0; i < num_qubits; i++){
    first_row_array[i] = -1;
  }
  bsp_push_reg(first_row_array, num_qubits * sizeof(int)); // Array to tell which processor has the lowest row. Used in swapping
  int *last_row_array = malloc(num_qubits * sizeof(int));
  for (int i = 0; i < num_qubits; i++){
    last_row_array[i] = -1;
  }
  bsp_push_reg(last_row_array, num_qubits * sizeof(int)); // Array to tell which processor has the highest row. Used in swapping

  bool empty = true;
  bsp_push_reg(&empty, sizeof(bool)); // Flag

  int *phase_summands = malloc(row_dim * M * sizeof(int));
  for (int i = 0; i < row_dim * M; i++){
    phase_summands[i] = -1;
  }
  bsp_push_reg(phase_summands, row_dim * M * sizeof(int));  // Array to calculate phase related calculations

  bsp_sync();

  // Code to obtain the first block of Hadamards. Extremely similar to Algorithm 1
  // ==================================================================================================================================================
  // ==================================================================================================================================================
  int k = 0;
  int diag_pid = 0;
  int row_w_one = 0;
  int last_row_w_one = 0;
  int temp_phase = 0;
  for (int col = num_qubits; col < 2 * num_qubits; col++){
    empty = true;
    diag_pid = k % N + (col % M) * N; // Processor that handles the right diagonal element

    // All processors that handle column "col" look to see if they have a 1 in that column
    if (col_name == col % M){
      for (int i = (int)floor(k / N); i < row_dim; i++){
        if (my_matrix[i][(int)floor(col / M)] == 1 && (i * N + row_name) >= k){
          empty = false;  // flag that an element exists
          row_w_one = i * N + row_name; // inverse mapping referring to the global mat
          bsp_put(diag_pid, &row_w_one, first_row_array, row_w_one * sizeof(int), sizeof(int));
          for (int ii = 0; ii < P; ii++){
            bsp_put(ii, &empty, &empty, 0, sizeof(bool)); // Update the value of "empty" for all processors
          }
          break;
        }
      }
    }
    bsp_sync();

    if (empty == false){
      // Processor handling the diagonal element decides on the processor which it will swap rows with (winner), and announces it to all other processors
      if (pid == diag_pid){
        for (int i = 0; i < num_qubits; i++){
          if (first_row_array[i] != -1){ // Any processor works. No need to play favorites 
            for (int ii = 0; ii < P; ii++){
              bsp_put(ii, &first_row_array[i], &winner_row, 0, sizeof(int));
            }
            break;
          }
        }
      }
      bsp_sync();

      // Swaps begin here! 
      // ============================================================================================================================================
      // Processors that control the row where the diagonal element currently is. They send the row to the processors responsible of the winner row
      if (row_name == k % N && row_name != winner_row % N){
        bsp_put(winner_row % N + (col_name % M) * N, my_matrix[(int)floor(k / N)], entries_received, 0, column_dim * sizeof(int));
        if (col_name == 0){
          bsp_put(winner_row % N + (col_name % M) * N, &my_vec[(int)floor(k / N)], &phase_received, 0, sizeof(int));
        }
      }

      //Processors that handle the winner row send to those handling row k.
      if (row_name == winner_row % N && row_name != k % N){
        bsp_put(k % N + (col_name % M) * N, my_matrix[(int)floor(winner_row / N)], entries_received, 0, column_dim * sizeof(int));
        if (col_name == 0){
          bsp_put(k % N + (col_name % M) * N, &my_vec[(int)floor(winner_row / N)], &phase_received, 0, sizeof(int));
        }
      }
      bsp_sync();

      // Processors overwrite their values of "my_matrix" with the elements just received
      // ------------------------------------------------------------------
      // This is the case when the processors handling row k also handle the winner row
      if (row_name == k % N && row_name == winner_row % N){
        int *temp_array = NULL;
        temp_array = my_matrix[(int)floor(k / N)];
        my_matrix[(int)floor(k / N)] = my_matrix[(int)floor(winner_row / N)];
        my_matrix[(int)floor(winner_row / N)] = temp_array;
        temp_array = NULL;

        // Swap phases
        if (col_name == 0){
          temp_phase = my_vec[(int)floor(k / N)];
          my_vec[(int)floor(k / N)] = my_vec[(int)floor(winner_row / N)];
          my_vec[(int)floor(winner_row / N)] = temp_phase;
          temp_phase = -1;
        }
      }
      // Case when the rocessors handling row k do not handle the winner row
      else if (row_name == k % N){ // Processors that handle row k update its local matrix
        for (int i = 0; i < column_dim; i++){
          my_matrix[(int)floor(k / N)][i] = entries_received[i];
        }
        if (col_name == 0){
          my_vec[(int)floor(k / N)] = phase_received;
        }
      }
      else if (row_name == winner_row % N){ // Processors that handle winner row update its local matrix
        for (int i = 0; i < column_dim; i++){
          my_matrix[(int)floor(winner_row / N)][i] = entries_received[i];
        }
        if (col_name == 0){
          my_vec[(int)floor(winner_row / N)] = phase_received;
        }
      }
    }
    // Search backwards on the left hand-side of the global matrix to obtain k2. And swap rows
    else{
      // Processors that have a Z literal (value 1 in "left_side_col", value 0 in "col") will set a flag and send it to "diag_pid"
      int left_side_col = col - num_qubits;
      if (col_name == left_side_col % M){
        empty = true;
        diag_pid = k % N + (left_side_col % M) * N;
        for (int i = num_qubits - 1; i >= k; i--){          // The loop starts from the back
          // All processors look to see if they have a Z literal
          if (my_matrix[(int)floor(i / N)][(int)floor(left_side_col / M)] == 1 && my_matrix[(int)floor(i / N)][(int)floor((left_side_col + num_qubits) / M)] == 0){
            if (((int)floor(i / N)) * N + row_name >= k){
              empty = false;
              last_row_w_one = ((int)floor(i / N)) * N + row_name;
              bsp_put(diag_pid, &last_row_w_one, last_row_array, last_row_w_one * sizeof(int), sizeof(int));
              for (int ii = 0; ii < P; ii++){
                bsp_put(ii, &empty, &empty, 0, sizeof(bool)); // Set the flag
              }
              break;
            }
          }
        }
      }
      bsp_sync();

      // If a Z literal was found
      if (empty == false){
        // Processor handling the diagonal element announces winner row
        if (pid == diag_pid){
          for (int i = num_qubits - 1; i >= k; i--){
            if (last_row_array[i] != -1){
              for (int ii = 0; ii < P; ii++){
                bsp_put(ii, &last_row_array[i], &winner_row, 0, sizeof(int));
              }
              break;
            }
          }
        }
        bsp_sync();

        // Processors that control the row where the diagonal element currently is. They send the row to the processors responsible of the winner row
        if (row_name == k % N && row_name != winner_row % N){
          bsp_put(winner_row % N + (col_name % M) * N, my_matrix[(int)floor(k / N)], entries_received, 0, column_dim * sizeof(int));
          if (col_name == 0){
            bsp_put(winner_row % N + (col_name % M) * N, &my_vec[(int)floor(k / N)], &phase_received, 0, sizeof(int));
          }
        }


        // Processors that handle the winner row send to those handling row k
        if (row_name == winner_row % N && row_name != k % N){
          bsp_put(k % N + (col_name % M) * N, my_matrix[(int)floor(winner_row / N)], entries_received, 0, column_dim * sizeof(int));
          if (col_name == 0){
            bsp_put(k % N + (col_name % M) * N, &my_vec[(int)floor(winner_row / N)], &phase_received, 0, sizeof(int));
          }
        }
        bsp_sync();

        // Processors overwrite their values of "my_matrix" with the elements just received
        // ------------------------------------------------------------------
        // This is the case when the processors handling row k also handle the winner row
        if (row_name == k % N && row_name == winner_row % N){
          // Swap rows
          int *temp_array = NULL;
          temp_array = my_matrix[(int)floor(k / N)];
          my_matrix[(int)floor(k / N)] = my_matrix[(int)floor(winner_row / N)];
          my_matrix[(int)floor(winner_row / N)] = temp_array;
          temp_array = NULL;

          // Swap phases
          if (col_name == 0){
            temp_phase = my_vec[(int)floor(k / N)];
            my_vec[(int)floor(k / N)] = my_vec[(int)floor(winner_row / N)];
            my_vec[(int)floor(winner_row / N)] = temp_phase;
            temp_phase = -1;
          }
        }
        // Case when the rocessors handling row k do not handle the winner row
        else if (row_name == k % N){ // Processors that handle row k update its local matrix
          for (int i = 0; i < column_dim; i++){
            my_matrix[(int)floor(k / N)][i] = entries_received[i];
          }
          if (col_name == 0){
            my_vec[(int)floor(k / N)] = phase_received;
          }
        }
        else if (row_name == winner_row % N){ // Processors that handle winner row update its local matrix
          for (int i = 0; i < column_dim; i++){
            my_matrix[(int)floor(winner_row / N)][i] = entries_received[i];
          }
          if (col_name == 0){
            my_vec[(int)floor(winner_row / N)] = phase_received;
          }
        }

        // Now we perform the conjugation! For the Hadamard this is just column swap
        // ----------------------------------------------------------------------------
        // First, processors in charge of row "k" mark where they have found an X,Y, or Z literal 
        if (row_name == k % N){
          for (int j = col + 1; j < 2 * num_qubits; j++){
            if (my_matrix[(int)floor(k / N)][(int)floor(j / M)] == 1 || my_matrix[(int)floor(k / N)][(int)floor((j - num_qubits) / M)] == 1){
              if (((int)floor(j / M)) * M + col_name >= j){
                has_one[((int)floor(j / M)) * M + col_name] = true; // Globally flag when there is a one in row k
              }
            }
          }
        }
        bsp_sync();

        // Create the circuit for the first block of Hadamards
        if (pid == 0){
          for (int i = num_qubits; i < 2 * num_qubits; i++){
            if (has_one[i] == 1){
              h1append(i - num_qubits);
            }
          }
        }

        int *loc_summands = malloc(row_dim * sizeof(int));
        for (int i = 0; i < row_dim; i++){
          loc_summands[i] = 0;
        }
        // All processors smartly look through "has_one" array and swap rows (locally!) if they need to
        for (int j = num_qubits + col_name; j < 2 * num_qubits; j = j + M){ 
          if (has_one[j] == 1){
            for (int i = 0; i < row_dim; i++){ // Loop through all rows, swap elements in row
              int temp_num = my_matrix[i][(int)floor((j - num_qubits) / M)];
              my_matrix[i][(int)floor((j - num_qubits) / M)] = my_matrix[i][(int)floor(j / M)];
              my_matrix[i][(int)floor(j / M)] = temp_num;

              loc_summands[i] = loc_summands[i] + my_matrix[i][(int)floor(j / M)] * my_matrix[i][(int)floor((j - num_qubits) / M)];
            }
          }
        }
        bsp_put(row_name, loc_summands, phase_summands, col_name * row_dim * sizeof(int), row_dim * sizeof(int));
        bsp_sync();

        // Update the phases
        if (col_name == 0){
          for (int i = 0; i < row_dim; i++){
            int temp_summands = 0;
            for (int l = i; l < row_dim * M; l = l + row_dim){
              temp_summands = temp_summands + phase_summands[l];
            }
            my_vec[i] = (my_vec[i] + temp_summands) % 2;
          }
        }
      }
    }

    // Reseting vars before the next iteration
    for (int i = 0; i < P; i++){
      winner_row_array[i] = -1;
    }
    for (int i = 0; i < column_dim; i++){
      entries_received[i] = -1;
    }
    for (int i = 0; i < 2 * num_qubits; i++){
      has_one[i] = false;
    }
    for (int i = 0; i < num_qubits; i++){
      last_row_array[i] = -1;
      first_row_array[i] = -1;
    }
    phase_received = -1;
    last_row_w_one = -1;
    row_w_one = -1;

    k++;
  }

  // Code to obtain the CNOT block
  // Instead of having communication between processors, we let each processor update the global matrix. Others can obtain the data necessary from there
  // =====================================================================================================================================
  // ====================================================================================================================================
  for (int row = 0; row < num_qubits; row++){ // We iterate throught the rows now!
    for (int k_ = row + num_qubits + 1; k_ < 2 * num_qubits; k_++){ // Look at right side of "mat" only
      if (row_name == row % N && my_matrix[(int)floor(row / N)][(int)floor(k_ / M)] == 1){ // if there are X or Y literals
        if (((int)floor(k_ / M)) * M + col_name >= k_){
          has_one[((int)floor(k_ / M)) * M + col_name] = true;
          found_a_one = true; // Set global flag
        }
      }
    }
    bsp_sync();

    if (found_a_one == true){
      // Processors read the has_one array and if 1, they update the columns in "mat" with their values
      for (int j = num_qubits + col_name; j < 2 * num_qubits; j = j + M){
        if (has_one[j] == 1){
          for (int i = 0; i < row_dim; i++){
            mat[i * N + row_name][j] = my_matrix[i][(int)floor(j / M)]; // Update the x part
            mat[i * N + row_name][j - num_qubits] = my_matrix[i][(int)floor((j - num_qubits) / M)]; // Update the z part
          }
        }
      }

      // Processors whose columns are indexed by "row" update its columns in mat (regardless of whether ones were found or not)
      if (col_name == row % M){
        for (int i = 0; i < row_dim; i++){
          mat[i * N + row_name][row + num_qubits] = my_matrix[i][(int)floor((row + num_qubits) / M)]; // Update the x part
          mat[i * N + row_name][row] = my_matrix[i][(int)floor(row / M)]; // Update the z part
        }
      }
      bsp_sync();

      // Update the phases
      if (col_name == 0){
        int x_row = 0;
        int z_row = 0;
        int x_j = 0;
        int z_j = 0;
        for (int j = num_qubits; j < 2 * num_qubits; j++){
          if (has_one[j] == 1){
            for (int l = 0; l < row_dim; l++){
              x_row = mat[l * N + row_name][row + num_qubits];
              z_row = mat[l * N + row_name][row];
              x_j = mat[l * N + row_name][j];
              z_j = mat[l * N + row_name][j - num_qubits];

              my_vec[l] = (my_vec[l] + x_row * z_j * (x_j + z_row + 1));
            }
          }
        }
      }

      // Processors whose columns are indexed by "row" update their local matrix. These have to search through the whole array "has_one"
      if (col_name == row % M){
        for (int j = num_qubits; j < 2 * num_qubits; j++){
          if (has_one[j] == 1){
            for (int i = 0; i < row_dim; i++){
              my_matrix[i][(int)floor(row / M)] = (my_matrix[i][(int)floor(row / M)] + mat[i * N + row_name][j - num_qubits]) % 2; // z_row += z_j
            }
          }
        }
      }

      // All processors smartly search through "has_one" and update when necessary
      for (int j = num_qubits + col_name; j < 2 * num_qubits; j = j + M){
        if (has_one[j] == 1){
          for (int i = 0; i < row_dim; i++){
            my_matrix[i][(int)floor(j / M)] = (my_matrix[i][(int)floor(j / M)] + mat[i * N + row_name][row + num_qubits]) % 2;
          }
        }
      }

      // Append to the circuit
      if (pid == 0){
        for (int j = num_qubits; j < 2 * num_qubits; j++){
          if (has_one[j] == 1){
            cxappend(row, j - num_qubits);
          }
        }
      }
      bsp_sync();

      // Reset vars
      if (pid == 0){
        for (int i = num_qubits; i < 2 * num_qubits; i++){
          has_one[i] = false;
        }
        found_a_one = false;
      }
    }
    bsp_sync();

  }

  // Code to obtain the CZ block
  // =====================================================================================================================================
  // =====================================================================================================================================
  for (int row = 0; row < num_qubits; row++){
    // Processors flag when they find a one in row row
    for (int k_ = (int)floor((row + 1) / M); k_ < column_dim / 2; k_++){ // Look at right side of mat
      if (row_name == row % N && my_matrix[(int)floor(row / N)][k_] == 1 && my_matrix[(int)floor(row / N)][k_ + column_dim / 2] == 0){ // If they are Z's
        if ((k_ * M + col_name) >= row + 1){
          has_one[(k_ + column_dim / 2) * M + col_name] = true; // Put on the right hand side of "has_one"
          found_a_one = true; // Set the global flag
        }
      }
    }
    bsp_sync();

    if (found_a_one == true){
      // Processors read "has_one" and if one, then they update the columns in "mat" with their values
      for (int j = column_dim / 2; j < column_dim; j++){
        if (has_one[j * M + col_name] == 1){
          for (int i = 0; i < row_dim; i++){
            mat[i * N + row_name][j * M + col_name] = my_matrix[i][j];  // Update the x part
            mat[i * N + row_name][(j * M + col_name) - num_qubits] = my_matrix[i][j - column_dim / 2]; // Update the z part
          }
        }
      }

      // Processors whose columns are indexed by "row", update its columns in mat (regardless of whether ones were found or not)
      if (col_name == row % M){
        for (int i = 0; i < row_dim; i++){
          mat[i * N + row_name][row + num_qubits] = my_matrix[i][(int)floor((row + num_qubits) / M)]; // Update the x part
          mat[i * N + row_name][row] = my_matrix[i][(int)floor(row / M)]; // Update the z part
        }
      }
      bsp_sync();

      // Processors whose columns are indexed by "row" update their local matrix. These have to search through the whole array "has_one"
      if (col_name == row % M){
        for (int j = num_qubits; j < 2 * num_qubits; j++){
          if (has_one[j] == 1){
            for (int i = 0; i < row_dim; i++){
              my_matrix[i][(int)floor(row / M)] = (my_matrix[i][(int)floor(row / M)] + mat[i * N + row_name][j]) % 2; // z_r += x_k'
            }
          }
        }
      }

      // All processors smartly search through "has_one" and update when necessary
      for (int j = num_qubits + col_name; j < 2 * num_qubits; j = j + M){
        if (has_one[j] == 1){
          for (int i = 0; i < row_dim; i++){
            my_matrix[i][(int)floor((j - num_qubits) / M)] = (my_matrix[i][(int)floor((j - num_qubits) / M)] + mat[i * N + row_name][row + num_qubits]) % 2; // z_k' += x_r
          }
        }
      }

      // Update the phases
      if (col_name == 0){
        int x_row = 0;
        int x_j = 0;
        for (int j = num_qubits; j < 2 * num_qubits; j++){
          if (has_one[j] == 1){
            for (int l = 0; l < row_dim; l++){
              x_row = mat[l * N + row_name][row + num_qubits];
              x_j = mat[l * N + row_name][j];

              my_vec[l] = (my_vec[l] + x_row * x_j) % 2;
            }
          }
        }
      }

      // Append to the circuit
      if (pid == 0){
        for (int j = num_qubits; j < 2 * num_qubits; j++){
          if (has_one[j] == 1){
            czappend(row, j - num_qubits);
          }
        }
      }
      bsp_sync();

      if (pid == 0){
        for (int i = num_qubits; i < 2 * num_qubits; i++){
          has_one[i] = false;
        }
        found_a_one = false;
      }
    }
    bsp_sync();
  }

  // Code to obtain the S block
  // =====================================================================================================================================
  // =====================================================================================================================================
  for (int j = num_qubits; j < 2 * num_qubits; j++){
    //Processors flag when they find a one in row j.
    if (row_name == j % N && col_name == j % M){
      if (my_matrix[(int)floor((j - num_qubits) / N)][(int)floor(j / M)] == 1 && my_matrix[(int)floor((j - num_qubits) / N)][(int)floor((j - num_qubits) / M)] == 1){ // If there is a Y literal
        has_one[j] = true; // Put on the right hand side of "has_one"
        found_a_one = true; // Set the global flag
        sappend(j - num_qubits);  // Append to the circuit
      }
    }
    bsp_sync();

    if (found_a_one == true){
      // Now we perform the conjugation. For an S gate this is just z_j += x_j
      int *loc_summands = malloc(row_dim * sizeof(int));
      for (int i = 0; i < row_dim; i++){
        loc_summands[i] = 0;
      }
      if (col_name == j % M){
        // First we do phase calculations, then we update the matrix!
        for (int i = 0; i < row_dim; i++){
          loc_summands[i] = my_matrix[i][(int)floor((j - num_qubits) / M)] * my_matrix[i][(int)floor(j / M)];
          my_matrix[i][(int)floor((j - num_qubits) / M)] = (my_matrix[i][(int)floor((j - num_qubits) / M)] + my_matrix[i][(int)floor(j / M)]) % 2;
        }
        bsp_put(row_name, loc_summands, phase_summands, 0, row_dim * sizeof(int)); // No matter the processor, it puts in the beginning of phase_summands. No interference with others guaranteed
      }
      bsp_sync();

      // Now we update the phases
      if (col_name == 0){
        for (int i = 0; i < row_dim; i++){
          my_vec[i] = (my_vec[i] + phase_summands[i]) % 2;
        }
      }
    }

    // Reset vars
    if (pid == 0){
      has_one[j] = false;
      found_a_one = false;
    }
    bsp_sync();
  }

  // Code to obtain the second block of Hadamards
  // ==================================================================================================================================================
  // ==================================================================================================================================================
  for (int j = num_qubits; j < 2 * num_qubits; j++){
    // Processors flag when they find a one in row j
    if (row_name == j % N && col_name == j % M){
      if (my_matrix[(int)floor((j - num_qubits) / N)][(int)floor(j / M)] == 1 && my_matrix[(int)floor((j - num_qubits) / N)][(int)floor((j - num_qubits) / M)] == 0){ // If there is a X literal
        has_one[j] = true; // Put on the right hand side of "has_one"
        found_a_one = true; // Set the global flag
        h2append(j - num_qubits); // Append to the circuit
      }
    }
    bsp_sync();

    if (found_a_one == true){
      // Now we do the conjugation, for a Hadamard this is just a local column swap
      int *loc_summands = malloc(row_dim * sizeof(int));
      for (int i = 0; i < row_dim; i++){
        loc_summands[i] = 0;
      }
      if (col_name == j % M){
        for (int i = 0; i < row_dim; i++){
          loc_summands[i] = my_matrix[i][(int)floor((j - num_qubits) / M)] * my_matrix[i][(int)floor(j / M)];
          int temp_num = my_matrix[i][(int)floor((j - num_qubits) / M)];
          my_matrix[i][(int)floor((j - num_qubits) / M)] = my_matrix[i][(int)floor(j / M)];
          my_matrix[i][(int)floor(j / M)] = temp_num;
        }
        bsp_put(row_name, loc_summands, phase_summands, 0, row_dim * sizeof(int));
      }
      bsp_sync();

      // Now we update the phases
      if (col_name == 0){
        for (int i = 0; i < row_dim; i++){
          my_vec[i] = (my_vec[i] + phase_summands[i]) % 2;
        }
      }
    }

    // Reset vars
    if (pid == 0){
      has_one[j] = false;
      found_a_one = false;
    }
    bsp_sync();
  }

  // Eliminate trailing Z literals
  // ==================================================================================================================================================
  // ==================================================================================================================================================
  for (int j = 0; j < num_qubits; j++){
    if (col_name == j % M){
      for (int i = (int)floor((j + 1) / N); i < row_dim; i++){
        //Loop through all of my rows.
        if (my_matrix[i][(int)floor(j / M)] == 1 && my_matrix[i][(int)floor((j + num_qubits) / M)] == 0){ // If there is a Z literal
          if (i * N + row_name > j){
            has_one[i * N + row_name] = true;
            found_a_one = true;
          }
        }
      }
    }
    bsp_sync();

    if (found_a_one == true){
      // Processors that handle row j, update their values in mat
      if (row_name == j % N){
        for (int g = 0; g < column_dim; g++){
          mat[j][g * M + col_name] = my_matrix[(int)floor(j / N)][g];
        }
        if (col_name == 0){
          vec[j] = my_vec[(int)floor(j / N)];
        }
      }
      bsp_sync();

      // All processors look in "has_one" to see if they have to update their own matrix
      int symplectic_prod = 0;
      int norm_prod = 0; // norm_prod = b_loc in Algorithm 1
      int v = 0;
      for (int l = row_name; l < num_qubits; l = l + N){
        if (has_one[l] == 1){
          symplectic_prod = 0;
          norm_prod = 0;
          v = 0;
          if (col_name == 0){ // Update the phases by computing the symplectic product, as in Algorithm 1
            for (int g = 0; g < num_qubits; g++){
              // symplectic_prod = z_l * x_j - x_l * z_j
              symplectic_prod += mat[l][g] * mat[j][g + num_qubits] - mat[l][g + num_qubits] * mat[j][g];
              v = ((my_matrix[(int)floor(l / N)][(int)floor(g / M)] + mat[j][g]) * (my_matrix[(int)floor(l / N)][(int)floor((g + num_qubits) / M)] + mat[j][g + num_qubits])) % 2;
              norm_prod += v - (my_matrix[(int)floor(l / N)][(int)floor(g / M)] + mat[j][g]) * (my_matrix[(int)floor(l / N)][(int)floor((g + num_qubits) / M)] + mat[j][g + num_qubits]);
            }
            my_vec[(int)floor(l / N)] = (MOD(symplectic_prod + norm_prod, 4) / 2 + my_vec[(int)floor(l / N)] + vec[j]) % 2;
          }
          for (int g = 0; g < column_dim; g++){
            my_matrix[(int)floor(l / N)][g] = (my_matrix[(int)floor(l / N)][g] + mat[j][g * M + col_name]) % 2;
          }
        }
      }
      bsp_sync();
    }

    // Reset vars
    if (pid == 0){
      found_a_one = false;
      for (int l = 0; l < num_qubits; l++){
        has_one[l] = false;
      }
    }
    bsp_sync();
  }

  // Delete the registers used for communication between processors
  bsp_pop_reg(&phase_received);
  bsp_pop_reg(entries_received);
  bsp_pop_reg(&empty);
  bsp_pop_reg(winner_row_array);
  bsp_pop_reg(last_row_array);
  bsp_pop_reg(&winner_row);
  bsp_pop_reg(phase_summands);
  bsp_pop_reg(first_row_array);

  // From their local matrix, all processors together reconstruct the global matrix and the global vector of phases
  for (int i = 0; i < row_dim; i++){
    for (int j = 0; j < column_dim; j++){
      mat[i * N + row_name][j * M + col_name] = my_matrix[i][j];
    }
    if (col_name == 0){
      vec[i * N + row_name] = my_vec[i];
    }
  }

  bsp_end();

  // Print the final matrix and circuit of gates
  for(int i = 0; i < num_qubits; i++){
    for(int j = 0; j < 2*num_qubits; j++){
      printf("%d,", mat[i][j]);
    }
    printf("\n");
  }
  for(int i =0; i < num_qubits; i++){
    printf("%d,", vec[i]);
  }
  for(int i =0; i < h1counter; i++){
    printf("(%s,%d,%d)", had1[i].name, had1[i].i, had1[i].j);
  }
  for(int i =0; i < cxcounter; i++){
    printf("(%s,%d,%d)", cnot[i].name, cnot[i].i, cnot[i].j);
  }
  for(int i =0; i < czcounter; i++){
    printf("(%s,%d,%d)", cphase[i].name, cphase[i].i, cphase[i].j);
  }
  for(int i =0; i < scounter; i++){
    printf("(%s,%d,%d)", phase[i].name, phase[i].i, phase[i].j);
  }
  for(int i =0; i < h2counter; i++){
    printf("(%s,%d,%d)", had2[i].name, had2[i].i, had2[i].j);
  }

  // Write the final matrix into a file
  write_output();
}

int main(int argc, char **argv){

  // Declare that "parallelalg1" is our parallel function
  bsp_init(parallelalg2, argc, argv);

  // Read a matrix from CSV and the vector of phases
  read_mat();
  read_vec();

  // Call the parallel function
  parallelalg2();
  
  exit(EXIT_SUCCESS);
}
