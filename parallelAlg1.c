/*Function that performs the first the canonical reduction of the stabilizer tableau (Algorithm 1 in the paper)

Constraints:
 - We have to make sure that N>=M
 - num_qubits % N = 0

 Inputs:
 - A valid stabilizer tableau. That is all of the rows are linearly independent, and contain at least a one
 - A vector of phases corresponding to the rows of the tableau. 

 Outputs:
 - A csv file containing the canonical decomposition of the tableau
*/ 

#include "stdio.h"
#include "stdlib.h"
#include "bsp.h"
#include "math.h"
#include "string.h"
#include "stdbool.h"

#define MOD(a, b) (((a % b) + b) % b)

const int P = 4;      // Number of processors
const int N = 2;          // Processor rows
const int M = 2;          // Processor columns
const int num_qubits = 8; // Sometimes referred to as n
int our_nan = -1 * ((2 * num_qubits) + 1);
int mat[8][2 * 8]; // Size of the tableau
int vec[8];        // Size of the phase vector

// Function that reads the starting tableau and loads it into directly into our 2-d array "mat"
void read_csv(){

  char buffer[100000];
  char *record, *line;
  int i = 0, j = 0;
  FILE *fstream = fopen("~/ParallelStabilizerInnerProduct/Inputs/starting_tableau.csv", "r");
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

  char buffer[1024];
  char *record, *line;
  int i = 0;
  FILE *fstream = fopen("~/ParallelStabilzerInnerProduct/Inputs/starting_phases.csv", "r");
  if (fstream == NULL){
    printf("\n file opening failed ");
  }
  while ((line = fgets(buffer, 1024, fstream)) != NULL){
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
  FILE *alg1_output_mat;
  FILE *alg1_output_vec;

  alg1_output_mat = fopen("~/ParallelStabilizerInnerProduct/Outputs/alg1_output_mat.txt", "w");

  if (alg1_output_mat != NULL){
    for (int i = 0; i < num_qubits; i++){
      for (int j = 0; j < 2 * num_qubits; j++){
        if (i == num_qubits - 1 && j == 2 * num_qubits - 1){
          fprintf(alg1_output_mat, "%d \n", mat[i][j]);
        }
        else{
          fprintf(alg1_output_mat, "%d,", mat[i][j]);
        }
      }
    }
  }
  else{
    printf("The file could not be opened");
  }

  fclose(alg1_output_mat);
  alg1_output_vec = fopen("~/ParallelStabilizerInnerProduct/Outputs/alg1_output_vec.txt", "w");

  if (alg1_output_vec != NULL){
    for (int l = 0; l < num_qubits; l++){
      if (l != num_qubits - 1){
        fprintf(alg1_output_vec, "%d,", vec[l]);
      }
      else{
        fprintf(alg1_output_vec, "%d", vec[l]);
      }
    }
  }

  fclose(alg1_output_vec);
}

// Parallel function
void parallelalg1(){

  // Begin using "P" processors
  bsp_begin(P);

  // Name the processors using 2d naming 
  int pid = bsp_pid(); // P0 = 00, P1 = 10 P2 = 01, P3 = 11
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
  // Initialize "total_as" and fill with zeros
  int *my_vec = NULL;
  int *total_as = NULL;
  if (col_name == 0){
    my_vec = malloc(row_dim * sizeof(int));
    total_as = malloc(row_dim * sizeof(int));
    for (int i = 0; i < row_dim; i++){
      my_vec[i] = vec[row_name + i * N];
      total_as[i] = 0;
    }
  }

  // Initialize registers that enable data communication between the processors
  int *entries_received = malloc(column_dim * sizeof(int));
  bsp_push_reg(entries_received, column_dim * sizeof(int)); // Array for received row entries. Used in swapping
  int phase_received = -1;
  bsp_push_reg(&phase_received, sizeof(int)); // Array for received phase entries
  int *winner_row_array = malloc(P * sizeof(int)); 
  for (int i = 0; i < P; i++){
    winner_row_array[i] = -1;
  }
  bsp_push_reg(winner_row_array, P * sizeof(int)); // Array to determine which processor wins
  int winner_row = -1;
  bsp_push_reg(&winner_row, sizeof(int)); // Which row wins
  bool empty = true;
  bsp_push_reg(&empty, sizeof(bool)); // Flag
  int *flag_array = malloc(row_dim * sizeof(int));
  for (int i = 0; i < row_dim; i++){
    flag_array[i] = -1;
  }
  bsp_push_reg(flag_array, row_dim * sizeof(int)); // Array for received column entries
  int *num_to_add = malloc(column_dim * sizeof(int));
  for (int i = 0; i < column_dim; i++){
    num_to_add[i] = -1;
  }
  bsp_push_reg(num_to_add, column_dim * sizeof(int)); // Array for received row entries. Used in row reduction

  int *partial_as = malloc(row_dim * M * sizeof(int));
  for (int i = 0; i < row_dim * M; i++){
    partial_as[i] = -1;
  }
  bsp_push_reg(partial_as, row_dim * M * sizeof(int)); // Helper array for phase calculations

  int *partial_bs = malloc(row_dim * M * sizeof(int));
  for (int i = 0; i < row_dim * M; i++){
    partial_bs[i] = -1;
  }
  bsp_push_reg(partial_bs, row_dim * M * sizeof(int)); // Helper array for phase calculations

  bsp_sync();

  // Helper local vars
  int *a_loc = malloc(row_dim * sizeof(int));
  int *b_loc = malloc(row_dim * sizeof(int));
  int k = 0;
  int diag_pid = 0;
  int row_w_one = 0;
  int temp_phase = 0;

  // BEGINNING OF THE MOST OUTER LOOP. REDUCTION FOR THE RIGHT SIDE OF THE TABLEAU
  // =====================================================================================================================================================
  // =====================================================================================================================================================
  // If there exists a one in a given column for all columns, then col=k. But if there isn't they will differ
  for (int col = num_qubits; col < 2 * num_qubits; col++){
    empty = true;

    diag_pid = k % N + (col % M) * N; // pid of the processor that handles the diagonal element (right side of the tableau)

    // All processors search for a 1 in column col. If they do have a 1, they communicate the number of the row to "diag_pid"
    if (col_name == col % M){
      for (int i = 0; i < row_dim; i++){
        if (my_matrix[i][(int)floor(col / M)] == 1 && (i * N + row_name) >= k){
          empty = false;
          row_w_one = i * N + row_name; // inverse mapping referring to the global mat
          bsp_put(diag_pid, &row_w_one, winner_row_array, pid * sizeof(int), sizeof(int));
          for (int ii = 0; ii < P; ii++){ // Update the value of "empty" for all processors
            bsp_put(ii, &empty, &empty, 0, sizeof(bool));
          }
          break;
        }
      }
    }
    bsp_sync();

    if (empty == false){
      // Processor handling the diagonal element decides on the processor which it will swap rows with (winner), and announces it to all other processors
      if (pid == diag_pid){
        if (winner_row_array[diag_pid] != -1){ // If "diag_pid" has the row with one, it should be the winner. Reduces communication costs
          for (int i = 0; i < P; i++){
            bsp_put(i, &winner_row_array[pid], &winner_row, 0, sizeof(int)); // Announce the winner row to all other processors
          }
        }
        else{
          for (int i = 0; i < P; i++){
            if (winner_row_array[i] != -1){
              for (int ii = 0; ii < P; ii++){
                bsp_put(ii, &winner_row_array[i], &winner_row, 0, sizeof(int)); // Announce the winner row to all other processors
              }
              break;
            }
          }
        }
      }
      bsp_sync();

      // Swaps begin here! 
      // ============================================================================================================================================
      // Processors that control the row where the diagonal element currently is. They send the row to the processors responsible of the winner row
      if (row_name == k % N && row_name != winner_row % N){
        bsp_put(winner_row % N + (col_name % M) * N, my_matrix[(int)floor(k / N)], entries_received, 0, column_dim * sizeof(int)); // Send row elements
        if (col_name == 0){
          bsp_put(winner_row % N + (col_name % M) * N, &my_vec[(int)floor(k / N)], &phase_received, 0, sizeof(int));  // Send phase elements
        }
      }

      // Processors that handle the winner row send to those handling row k
      if (row_name == winner_row % N && row_name != k % N){
        bsp_put(k % N + (col_name % M) * N, my_matrix[(int)floor(winner_row / N)], entries_received, 0, column_dim * sizeof(int)); // Send row elements
        if (col_name == 0){
          bsp_put(k % N + (col_name % M) * N, &my_vec[(int)floor(winner_row / N)], &phase_received, 0, sizeof(int));  // Send phase elements
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
      // Case when the processors handling row k do not handle the winner row
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
      bsp_sync();

      // Row reduction steps! Get rid of ones in the same column as "col"
      // =========================================================================================================================================
      // Processors that control row "k" send their row elements to the processor column
      if (row_name == k % N){
        for (int l = 0; l < N; l++){
          bsp_put(l + (col_name % M) * N, my_matrix[(int)floor(k / N)], num_to_add, 0, column_dim * sizeof(int));
          if (col_name == 0){
            // Processor in row k with column name 0 send their phases to all other procs with col_name 0. We reuse the "phase_received" memory spaces
            bsp_put(l, &my_vec[(int)floor(k / N)], &phase_received, 0, sizeof(int));
          }
        }
      }
      
      // Processors that control column "col" send their column entries to their processor row. We need two cases!
      // Case 1: processor handling column "col" and row "k" has to be careful, we do not want to delete its 1
      if (col_name == col % M && row_name == k % N){
        my_matrix[(int)floor(k / N)][(int)floor(col / M)] = -1;
        for (int i = 0; i < row_dim; i++){
          for (int l = 0; l < M; l++){
            bsp_put(row_name + (l % M) * N, &my_matrix[i][(int)floor(col / M)], flag_array, i * sizeof(int), sizeof(int));
          }
        }
        my_matrix[(int)floor(k / N)][(int)floor(col / M)] = 1;
      }
      // Case 2: everything proceeds normally
      else if (col_name == col % M){
        for (int i = 0; i < row_dim; i++){
          for (int l = 0; l < M; l++){
            bsp_put(row_name + (l % M) * N, &my_matrix[i][(int)floor(col / M)], flag_array, i * sizeof(int), sizeof(int));
          }
        }
      }
      bsp_sync();

      // Compute the first part of the symplectic inner product in parallel. All local!
      // Let z,x = my_matrix[:column_dim/2], my_matrix[column_dim/2:] and z',x' = nums_to_add[:column_dim/2], nums_to_add[column_dim/2:]
      // Compute a_loc = <z,x'> - <x,z'> and b_loc = <z+z', x+x'> % 2 - <z+z', x+x'>
      for (int l = 0; l < row_dim; l++){
        a_loc[l] = our_nan;
        b_loc[l] = our_nan;
      }
      for (int i = 0; i < row_dim; i++){
        if (flag_array[i] == 1){
          a_loc[i] = 0;
          b_loc[i] = 0;
          for (int j = 0; j < (int)(column_dim / 2); j++){
            a_loc[i] += (my_matrix[i][j] * num_to_add[j + (column_dim / 2)]) - (my_matrix[i][j + (column_dim / 2)] * num_to_add[j]);
            int v = MOD((my_matrix[i][j] + num_to_add[j]) * (my_matrix[i][j + column_dim / 2] + num_to_add[j + column_dim / 2]), 2);
            b_loc[i] += v - ((my_matrix[i][j] + num_to_add[j]) * (my_matrix[i][j + column_dim / 2] + num_to_add[j + column_dim / 2]));
          }
        }
      }
      // We put the local "a_loc" and "b_loc", on the processor with column_name = 0 that also handles that rows
      bsp_put(row_name, a_loc, partial_as, col_name * row_dim * sizeof(int), row_dim * sizeof(int));
      bsp_put(row_name, b_loc, partial_bs, col_name * row_dim * sizeof(int), row_dim * sizeof(int));
      bsp_sync();

      // The processor with column name 0, will use the partial a's to compute a for every row it particates in
      if (col_name == 0){
        for (int i = 0; i < row_dim; i++){
          if (partial_as[i] != our_nan){
            for (int ii = i; ii < row_dim * M; ii += row_dim){
              total_as[i] += partial_as[ii] + partial_bs[ii];
            }
            total_as[i] = MOD(total_as[i], 4);
            my_vec[i] = MOD(total_as[i] / 2 + my_vec[i] + phase_received, 2);
          }
        }
      }

      // Each processor replaces its rows and phases with the things just computed
      for (int i = 0; i < num_qubits; i++){
        if (row_name == i % N){
          if (flag_array[(int)floor(i / N)] == 1){
            for (int l = 0; l < column_dim; l++){
              my_matrix[(int)floor(i / N)][l] = (my_matrix[(int)floor(i / N)][l] + num_to_add[l]) % 2;
            }
          }
        }
      }

      // Reset vars before next iteration
      for (int i = 0; i < P; i++){
        winner_row_array[i] = -1;
      }
      for (int i = 0; i < row_dim; i++){
        flag_array[i] = -1;
        if (col_name == 0){
          total_as[i] = 0;
        }
      }
      for (int i = 0; i < column_dim; i++){
        num_to_add[i] = -1;
        entries_received[i] = -1;
      }
      for (int i = 0; i < row_dim * M; i++){
        partial_as[i] = our_nan;
        partial_bs[i] = our_nan;
      }
      phase_received = -1;

      k++;

      bsp_sync();
    }
  }

  // BEGINNING OF THE MOST OUTER LOOP. REDUCTION FOR THE LEFT SIDE OF THE TABLEAU
  // =====================================================================================================================================================
  // =====================================================================================================================================================
  // This part follows the exact structure as in the right case, the only difference is that "k" starts from where it left off in the previous section.
  for (int col = 0; col < num_qubits; col++){
    empty = true;

    diag_pid = k % N + (col % M) * N;

    if (col_name == col % M){
      for (int i = 0; i < row_dim; i++){
        if (my_matrix[i][(int)floor(col / M)] == 1 && (i * N + row_name) >= k){
          empty = false;
          row_w_one = i * N + row_name;
          bsp_put(diag_pid, &row_w_one, winner_row_array, pid * sizeof(int), sizeof(int));
          for (int ii = 0; ii < P; ii++){
            bsp_put(ii, &empty, &empty, 0, sizeof(bool));
          }
          break;
        }
      }
    }
    bsp_sync();

    if (empty == false){
      if (pid == diag_pid){
        if (winner_row_array[pid] != -1){
          for (int i = 0; i < P; i++){
            bsp_put(i, &winner_row_array[pid], &winner_row, 0, sizeof(int));
          }
        }
        else{
          for (int i = 0; i < P; i++){
            if (winner_row_array[i] != -1){
              for (int ii = 0; ii < P; ii++){
                bsp_put(ii, &winner_row_array[i], &winner_row, 0, sizeof(int));
              }
              break;
            }
          }
        }
      }
      bsp_sync();

      // Swaps begin here! 
      // ============================================================================================================================================
      // Processors that control the row where the diagonal element we are checking atm is. We send the row to the processors responsible of the winner row.
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

      //We overwrite our values in "my_matrix" with those received.
      // This is the case where the processor handling row k also handles the winner row.
      if (row_name == k % N && row_name == winner_row % N){
        int *temp_array = NULL;
        temp_array = my_matrix[(int)floor(k / N)];
        my_matrix[(int)floor(k / N)] = my_matrix[(int)floor(winner_row / N)];
        my_matrix[(int)floor(winner_row / N)] = temp_array;
        temp_array = NULL;

        if (col_name == 0){
          temp_phase = my_vec[(int)floor(k / N)];
          my_vec[(int)floor(k / N)] = my_vec[(int)floor(winner_row / N)];
          my_vec[(int)floor(winner_row / N)] = temp_phase;
          temp_phase = -1;
        }
      }
      else if (row_name == k % N){ //Processor that handles row k updates its matrix.
        for (int i = 0; i < column_dim; i++){
          my_matrix[(int)floor(k / N)][i] = entries_received[i];
        }
        if (col_name == 0){
          my_vec[(int)floor(k / N)] = phase_received;
        }
      }
      else if (row_name == winner_row % N){ //Processor that handles winner row updates its matrix.
        for (int i = 0; i < column_dim; i++){
          my_matrix[(int)floor(winner_row / N)][i] = entries_received[i];
        }
        if (col_name == 0){
          my_vec[(int)floor(winner_row / N)] = phase_received;
        }
      }

      // Row reduction steps! Get rid of ones in the same column as "col"
      // =========================================================================================================================================
      if (row_name == k % N){
        for (int l = 0; l < N; l++){
          bsp_put(l + (col_name % M) * N, my_matrix[(int)floor(k / N)], num_to_add, 0, (int)(column_dim / 2) * sizeof(int));
          if (col_name == 0){
            bsp_put(l, &my_vec[(int)floor(k / N)], &phase_received, 0, sizeof(int));
          }
        }
      }

      if (col_name == col % M && row_name == k % N){
        my_matrix[(int)floor(k / N)][(int)floor(col / M)] = -1;
        for (int i = 0; i < row_dim; i++){
          for (int l = 0; l < M; l++){
            bsp_put(row_name + (l % M) * N, &my_matrix[i][(int)floor(col / M)], flag_array, i * sizeof(int), sizeof(int));
          }
        }
        my_matrix[(int)floor(k / N)][(int)floor(col / M)] = 1;
      }
      else if (col_name == col % M){
        for (int i = 0; i < row_dim; i++){
          for (int l = 0; l < M; l++){
            bsp_put(row_name + (l % M) * N, &my_matrix[i][(int)floor(col / M)], flag_array, i * sizeof(int), sizeof(int));
          }
        }
      }
      bsp_sync();

      for (int l = 0; l < row_dim; l++){
        a_loc[l] = our_nan;
        b_loc[l] = our_nan;
      }
      for (int i = 0; i < row_dim; i++){
        if (flag_array[i] == 1){
          a_loc[i] = 0;
          b_loc[i] = 0;
          for (int j = 0; j < (int)column_dim / 2; j++){
            a_loc[i] += -1 * (my_matrix[i][j + (column_dim / 2)] * num_to_add[j]);
            int v = MOD((my_matrix[i][j] + num_to_add[j]) * my_matrix[i][j + column_dim / 2], 2);
            b_loc[i] += v - ((my_matrix[i][j] + num_to_add[j]) * my_matrix[i][j + column_dim / 2]);
          }
        }
      }
      bsp_put(row_name, a_loc, partial_as, col_name * row_dim * sizeof(int), row_dim * sizeof(int));
      bsp_put(row_name, b_loc, partial_bs, col_name * row_dim * sizeof(int), row_dim * sizeof(int));
      bsp_sync();

      if (col_name == 0){
        for (int i = 0; i < row_dim; i++){
          if (partial_as[i] != our_nan){
            for (int ii = i; ii < row_dim * M; ii += row_dim){
              total_as[i] += partial_as[ii] + partial_bs[ii];
            }
            total_as[i] = MOD(total_as[i], 4);
            my_vec[i] = MOD(total_as[i] / 2 + my_vec[i] + phase_received, 2);
          }
        }
      }

      for (int i = 0; i < num_qubits; i++){
        if (row_name == i % N){
          if (flag_array[(int)floor(i / N)] == 1){
            for (int l = 0; l < column_dim / 2; l++){
              my_matrix[(int)floor(i / N)][l] = (my_matrix[(int)floor(i / N)][l] + num_to_add[l]) % 2;
            }
          }
        }
      }

      for (int i = 0; i < P; i++){
        winner_row_array[i] = -1;
      }
      for (int i = 0; i < row_dim; i++){
        flag_array[i] = -1;
        if (col_name == 0){
          total_as[i] = 0;
        }
      }
      for (int i = 0; i < column_dim; i++){
        num_to_add[i] = -1;
        entries_received[i] = -1;
      }
      for (int i = 0; i < row_dim * M; i++){
        partial_as[i] = our_nan;
        partial_bs[i] = our_nan;
      }
      phase_received = -1;

      k++;

      bsp_sync();
    }
  }

  // Delete the registers used for communication between processors
  bsp_pop_reg(partial_as);
  bsp_pop_reg(partial_bs);
  bsp_pop_reg(&phase_received);
  bsp_pop_reg(entries_received);
  bsp_pop_reg(&empty);
  bsp_pop_reg(winner_row_array);
  bsp_pop_reg(&winner_row);
  bsp_pop_reg(flag_array);
  bsp_pop_reg(num_to_add);

  // From their local matrix, all processors together reconstruct the global matrix and the global vector of phases
  for (int j = 0; j < column_dim; j++){
      mat[i * N + row_name][j * M + col_name] = my_matrix[i][j];
  }
  if (col_name == 0){
      vec[i * N + row_name] = my_vec[i];
  }

  bsp_end();

  // Print the final matrix
  for (int i = 0; i < num_qubits; i++){
    for (int j = 0; j < 2 * num_qubits; j++){
      printf("%d,", mat[i][j]);
    }
    printf("\n");
  }

  // Write the final matrix into a file
  write_output();
}

int main(int argc, char **argv){
  
  // Declare that "parallelalg1" is our parallel function
  bsp_init(parallelalg1, argc, argv);

  // Read the starting stabilizer tableau and the vector of phases
  read_csv();
  read_vec();

  // Call the parallel function
  parallelalg1();

  exit(EXIT_SUCCESS);
}
