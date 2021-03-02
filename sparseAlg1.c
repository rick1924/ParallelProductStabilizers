/* Sparse function that performs the first the canonical reduction of the stabilizer tableau (Algorithm 1 in the paper). The non-zero elements are stored in a disordered doubly-linked list.

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
#include "string.h"
#include "math.h"
#include "bsp.h"
#include "stdbool.h"

const int P = 4;    // Number of processors
const int N = 2;    // Processor rows
const int M = 2;    // Processor columns
const int num_qubits = 8; // Sometimes referred to as n
int mat[8][2 * 8];  // Size of the tableau
int num_of_ones = 0;    // Counter for the number of ones in the tableau

// HEAD is the most important node of the list, as it is through this node that we can access the child nodes.
// Therefore, we must never delete it or change its address. If we loose it, then we cannot access the child
// nodes located in heap memory. HEAD is the link that connects different memory locations and keeps the list
// together.

// We define the doubly-linked list. It has two integer fields "row" and "column"
// "prev" and "next" point to the previous and next nodes respectively
struct node{
    int row;
    int column;
    struct node *prev;
    struct node *next;
};

// Function that creates nodes depending on several scenarios
// Note: if there is no relevant node in the list, HEAD is a fake node with entries "-1"
void append_to_list(struct node *head, int i, int j){
    // If HEAD is a fake node (list is empty)
    if (head->row == -1){
        head->row = i;
        head->column = j;
        head->next = NULL;
        head->prev = NULL;
    }
    // If there is at least one real node in list
    else{
        while (head->next != NULL){
            head = head->next;
        }
        // Create child node
        struct node *new = malloc(sizeof(struct node));
        new->row = i;
        new->column = j;
        head->next = new;
        new->prev = head;
        new->next = NULL;
    }
}

// Function that deletes a normal (non-HEAD) node
void delete_node(struct node *node){
    // If you are the only node
    if (node->next == NULL && node->prev == NULL){
        free(node);
    }
    // If the list ends with you
    else if (node->next == NULL){
        struct node *before = node->prev;
        before->next = NULL;
        free(node);
    }
    // If none of the above
    else{
        struct node *before = node->prev;
        struct node *after = node->next;
        before->next = node->next;
        after->prev = node->prev;
        free(node);
    }
}

// Function that "deletes" the HEAD node
void delete_node_head(struct node *head){
    // If there is no other node, then become a fake node
    if (head->next == NULL){
        head->row = -1;
        head->column = -1;
        head->prev = NULL;
    }
    // HEAD copies the information of the next node and we simply delete the next node
    else{
        struct node *after = head->next;
        head->next = after->next;
        head->row = after->row;
        head->column = after->column;
        if (after->next == NULL){
            free(after);
        }
        else{
            struct node *afterafter = after->next;
            afterafter->prev = head;
            free(after);
        }
    }
}

// Function that deletes any node.
void delete_single(struct node *head, int i, int j){
    struct node *temp = head;
    while (temp != NULL){
        if (temp->row == i && temp->column == j){
            if (temp == head){
                delete_node_head(temp); // Special function that "deletes" HEAD
                break;
            }
            else{
                delete_node(temp); // Delete a normal node
                break;
            }
        }
        temp = temp->next;
    }
}

// Function that deletes nodes that have the same two values, only keeps one node with these two values
void delete_repeated_pairs(struct node *head, int i, int j){
    struct node *temp = head;
    int count = 0;
    bool flag = false;
    while (temp != NULL){
        flag = false;
        if (temp->row != -1){
            if (temp->row == i && temp->column == j){
                if (temp == head){
                    delete_node_head(temp);
                    flag = true;
                }
                else{
                    delete_single(head, i, j);
                    flag = false;
                }
                count++;
            }
            if (flag == false){
                temp = temp->next;
            }
        }
        else{
            break;
        }
    }
    if (count == 1){
        append_to_list(head, i, j);
    }
}

// Function that deletes all nodes with a certain row value
void delete_multiple_row(struct node *head, int i){
    struct node *temp = head;
    struct node *temptemp;
    while (temp != NULL){
        if (temp->row == i){
            temptemp = temp;
            if (temptemp == head){
                delete_node_head(temptemp);
            }
            else{
                temp = temp->next;
                delete_node(temptemp);
            }
        }
        else{
            temp = temp->next;
        }
    }
}

// Function that deletes all nodes with a certain column value
void delete_multiple_column(struct node *head, int j){
    struct node *temp = head;
    struct node *temptemp;
    while (temp != NULL){
        if (temp->column == j){
            temptemp = temp;
            if (temptemp == head){
                delete_node_head(temptemp);
            }
            else{
                temp = temp->next;
                delete_node(temptemp);
            }
        }
        else{
            temp = temp->next;
        }
    }
}

// Find an a row that contains a one in column j, and that satisfies row >= roof
int find_i_given_j(struct node *head, int j, int roof){
    int found = -1;
    while (head != NULL){
        if (head->column == j && head->row >= roof){
            found = head->row;
            break;
        }
        head = head->next;
    }
    return found;
}

// Function that flags in which column there is a non-zero element, given a row i
int find_ones_row(struct node *head, int *array, int i){
    int counter = 0;
    while (head != NULL){
        if (head->row == i){
            array[counter] = head->column;
            counter++;
        }
        head = head->next;
    }
    return counter;
}

// Function that flags in which row there is a non-zero element, given a column j
int find_ones_col_except(struct node *head, int *array, int j, int num_except){
    int counter = 0;
    while (head != NULL){
        if (head->column == j && head->row != num_except){
            array[counter] = head->row;
            counter++;
        }
        head = head->next;
    }
    return counter;
}

// Function that swaps two rows
void swap_two_rows(struct node *head, int row1, int row2){
    while (head != NULL){
        if (head->row == row1){
            head->row = row2;
        }
        else if (head->row == row2){
            head->row = row1;
        }
        head = head->next;
    }
}

// Function that adds many rows together
// "flag_array" contains the rows that participate in the addition
// "num_to_add" contains the column values we should add to the list
// "sizeflag" and "sizenum" tell the size of the arrays
// Both "flag_array" and "num_to_add" may contain elements "-1" that tell the processor when to stop. Save time
void add_rows(struct node *head, int *flag_array, int *num_to_add, int sizeflag, int sizenum){
    bool flag = false;
    struct node *temp = head; // Set a temporary address
    for (int i = 0; i < sizeflag; i++){
        if (flag_array[i] != -1){ // Check if we can stop
            for (int j = 0; j < sizenum; j++){
                if (num_to_add[j] != -1){ // Check if we have to stop
                    flag = false;
                    while (temp != NULL){
                        if (temp->row == flag_array[i] && temp->column == num_to_add[j]){ // If element already exists, we delete it
                            delete_single(head, flag_array[i], num_to_add[j]);
                            flag = true;
                            break;
                        }
                        temp = temp->next;
                    }
                    // If element does not exist in the list we add it
                    if (flag == false){
                        append_to_list(head, flag_array[i], num_to_add[j]);
                    }
                    temp = head;
                }
                else{
                    break;
                }
            }
            temp = head;
        }
        else{
            break;
        }
    }
}

// Function that counts the number of nodes in the list
int count_tot_ones(struct node *head){
    int total = 0;
    while (head != NULL){
        total++;
        head = head->next;
    }
    return total;
}

// Function that prints the two relevant elements in the nodes
void print_list(struct node *head){
    if (head->next == NULL){
        printf("(%d, %d) \n", head->row, head->column);
    }
    else{
        while (head != NULL){
            printf("(%d, %d)", head->row, head->column);
            head = head->next;
        }
        printf("\n");
    }
}

// Function that reconstructs the global matrix from the linked-list
void reconstruct_mat(struct node *head){
    while (head != NULL){
        mat[head->row][head->column] = 1;
        head = head->next;
    }
}

// Function that reads the tableau in its canonical form and loads it into directly into our 2d array "mat"
void read_csv(){
    // char buffer[10000];
    char buffer[(num_qubits * num_qubits) * sizeof(int) + 1];
    char *record, *line;
    int i = 0, j = 0;
    FILE *fstream = fopen("~/ParallelStabilizerInnerProduct/Inputs/starting_tableau.csv", "r");
    // FILE *fstream = fopen("~/ParallelStabilizerInnerProduct/Inputs/sparse_tableau.csv", "r");
    if (fstream == NULL){
        printf("\n file opening failed ");
    }
    while ((line = fgets(buffer, sizeof(buffer), fstream)) != NULL){
        record = strtok(line, ",");
        while (record != NULL){
            mat[i][j++] = atoi(record);
            record = strtok(NULL, ",");
        }
        ++i;
    }

    fclose(fstream);
}

// Function that writes the tableau and resulting phases into a CSV file
void write_output(){
    FILE *sparsealg1_output_mat;

    sparsealg1_output_mat = fopen("~/ParallelStabilizerInnerProduct/Outputs/sparsealg1_output_mat.txt", "w");

    if (sparsealg1_output_mat != NULL){
        for (int i = 0; i < num_qubits; i++){
            for (int j = 0; j < 2 * num_qubits; j++){
                if (i == num_qubits - 1 && j == 2 * num_qubits - 1){
                    fprintf(sparsealg1_output_mat, "%d \n", mat[i][j]);
                }
                else{
                    fprintf(sparsealg1_output_mat, "%d,", mat[i][j]);
                }
            }
        }
    }
    else{
        printf("The file could not be opened");
    }

    fclose(sparsealg1_output_mat);
}

// Parallel function
void sparsealg1(){

    bsp_begin(P);   // Begin using "P" processors

    // Initialization
    // ------------------------------------------------------------------------
    //All processors make a fake node called head. This node is local
    int tot_ones = 0;
    struct node head;
    head.row = -1;
    head.column = -1;
    head.prev = NULL;
    head.next = NULL;

    // Name the processors using 2d naming 
    int pid = bsp_pid();    // P0 = 00, P1 = 10, P2 = 01, P3 = 11
    int row_name = pid % N;
    int col_name = (int)floor(pid / N);

    const int column_dim = 2 * num_qubits / M;
    const int row_dim = num_qubits / N;
    const int tot_values = column_dim * row_dim;

    // Each processor populates its doubly linked list
    for (int i = 0; i < row_dim; i++){
        for (int j = 0; j < column_dim; j++){
            if (mat[row_name + i * N][col_name + j * M] == 1){
                append_to_list(&head, i * N + row_name, j * M + col_name);
                tot_ones++;
            }
            mat[row_name + i * N][col_name + j * M] = 0;    // Turn the global matrix to all zeros
        }
    }

    // Update the global number of ones
    num_of_ones += tot_ones;

    // Initialize registers that enable data communication between the processors
    bool empty = true;
    bsp_push_reg(&empty, sizeof(bool));     // Flag
    int *winner_row_array = malloc(P * sizeof(int));
    for (int i = 0; i < P; i++){
        winner_row_array[i] = -1;
    }
    bsp_push_reg(winner_row_array, P * sizeof(int));    // Array to determine which processor wins
    int winner_row = -1;
    bsp_push_reg(&winner_row, sizeof(int));     // Which row wins
    int *entries_received = malloc(column_dim * sizeof(int));
    for (int i = 0; i < column_dim; i++){
        entries_received[i] = -1;
    }
    bsp_push_reg(entries_received, column_dim * sizeof(int));

    int *flag_array = malloc(row_dim * sizeof(int)); // Array that will tell the processors when they need to add or not.
    for (int i = 0; i < row_dim; i++){
        flag_array[i] = -1;
    }
    bsp_push_reg(flag_array, row_dim * sizeof(int));    // Array for received column entries
    int *num_to_add = malloc(column_dim * sizeof(int));
    for (int i = 0; i < column_dim; i++){
        num_to_add[i] = -1;
    }
    bsp_push_reg(num_to_add, column_dim * sizeof(int));     // Array for received row entries. Used in row reduction
    bsp_sync();

    int k = 0;
    int diag_pid = 0;
    int row_w_one = -1;
    int *all_ones_row = malloc(column_dim * sizeof(int)); // Array bounded by the number of columns they control
    int tot_ones_row = 0;
    int *all_ones_col = malloc(row_dim * sizeof(int));
    int tot_ones_col = 0;

  // BEGINNING OF THE MOST OUTER LOOP. REDUCTION FOR THE RIGHT SIDE OF THE TABLEAU
  // =====================================================================================================================================================
  // =====================================================================================================================================================
    for (int col = num_qubits; col < 2 * num_qubits; col++){

        // Loop dependent initializations
        for (int i = 0; i < column_dim; i++){
            all_ones_row[i] = -1;
        }
        for (int i = 0; i < row_dim; i++){
            all_ones_col[i] = -1;
        }
        empty = true;
        diag_pid = k % N + (col % M) * N;   // pid of the processor that handles the diagonal element (right side of the tableau)

        // All processor look through their list see if they have a one in column col
        if (col_name == col % M){
            row_w_one = find_i_given_j(&head, col, k);
            if (row_w_one != -1){
                empty = false;
                bsp_put(diag_pid, &row_w_one, winner_row_array, pid * sizeof(int), sizeof(int));
                for (int ii = 0; ii < P; ii++){     // Update the value of "empty" for all procs
                    bsp_put(ii, &empty, &empty, 0, sizeof(bool));
                }
            }
        }
        bsp_sync();

        if (empty == false){
            // Processor handling the diagonal element decides on the processor which it will swap rows with (winner), and announces it to all other processors
            if (pid == diag_pid){
                if (winner_row_array[diag_pid] != -1){ // If "diag_pid" has the row with one, it should be the winner. Reduces communication costs
                    for (int i = 0; i < P; i++){
                        bsp_put(i, &winner_row_array[pid], &winner_row, 0, sizeof(int));    // Announce the winner row to all other processors
                    }
                }
                else{
                    for (int i = 0; i < P; i++){
                        if (winner_row_array[i] != -1){
                            for (int ii = 0; ii < P; ii++){
                                bsp_put(ii, &winner_row_array[i], &winner_row, 0, sizeof(int));     // Announce the winner row to all other processors
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
                // These processors scan the list to find the ones in row "k". They put on the "winner_row" proc with the same column as them
                tot_ones_row = find_ones_row(&head, all_ones_row, k);
                bsp_put(winner_row % N + (col_name % M) * N, all_ones_row, entries_received, 0, tot_ones_row * sizeof(int));
                // Delete all of your ones in that row
                delete_multiple_row(&head, k);
            }

            // Processors that handle the winner row send to those handling row "k"
            if (row_name == winner_row % N && row_name != k % N){
                tot_ones_row = find_ones_row(&head, all_ones_row, winner_row);
                bsp_put(k % N + (col_name % M) * N, all_ones_row, entries_received, 0, tot_ones_row * sizeof(int));
                //Delete all of your ones in that row
                delete_multiple_row(&head, winner_row);
            }
            bsp_sync();

            // We update our list of nodes for the possible cases
            // In the case where the processor handles both rows, it only swaps the row entries of the list
            if (row_name == k % N && row_name == winner_row % N){
                swap_two_rows(&head, k, winner_row);
            }
            else if (row_name == k % N){
                for (int i = 0; i < column_dim; i++){
                    if (entries_received[i] != -1){
                        append_to_list(&head, k, entries_received[i]);
                    }
                    else{
                        break;
                    }
                }
            }
            else if (row_name == winner_row % N){
                for (int i = 0; i < column_dim; i++){
                    if (entries_received[i] != -1){
                        append_to_list(&head, winner_row, entries_received[i]);
                    }
                    else{
                        break;
                    }
                }
            }

            // Reset vars
            for (int i = 0; i < column_dim; i++){
                all_ones_row[i] = -1;
            }

            // Row reduction steps! Get rid of ones in the same column as "col"
            // =========================================================================================================================================
            // Processors that control row "k" send their row elements to the processor column
            if (row_name == k % N){
                tot_ones_row = find_ones_row(&head, all_ones_row, k);
                for (int l = 0; l < N; l++){
                    bsp_put(l + (col_name % M) * N, all_ones_row, num_to_add, 0, tot_ones_row * sizeof(int));
                }
            }

            // Processors that control column "col" send their column entries to their processor row. No need for two cases as in dense implementation
            if (col_name == col % M){
                tot_ones_col = find_ones_col_except(&head, all_ones_col, col, k);
                //Send to the processor row
                for (int l = 0; l < M; l++){
                    bsp_put(row_name + (l % M) * N, all_ones_col, flag_array, 0, tot_ones_col * sizeof(int));
                }
            }
            bsp_sync();

            // Do the updates
            add_rows(&head, flag_array, num_to_add, row_dim, column_dim);

            // Reset vars
            for (int i = 0; i < P; i++){
                winner_row_array[i] = -1;
            }
            for (int i = 0; i < column_dim; i++){
                entries_received[i] = -1;
                num_to_add[i] = -1;
            }
            for (int i = 0; i < row_dim; i++){
                flag_array[i] = -1;
            }

            k++;

            bsp_sync();
        }
    }

    // BEGINNING OF THE MOST OUTER LOOP. REDUCTION FOR THE LEFT SIDE OF THE TABLEAU
    // =====================================================================================================================================================
    // =====================================================================================================================================================
    for (int col = 0; col < num_qubits; col++){

        // Loop dependent initializations
        for (int i = 0; i < column_dim; i++){
            all_ones_row[i] = -1;
        }
        for (int i = 0; i < row_dim; i++){
            all_ones_col[i] = -1;
        }
        empty = true;
        diag_pid = k % N + (col % M) * N;

        // All processors look through their list see if they have a one in column "col"
        if (col_name == col % M){
            row_w_one = find_i_given_j(&head, col, k);
            if (row_w_one != -1){
                empty = false;
                bsp_put(diag_pid, &row_w_one, winner_row_array, pid * sizeof(int), sizeof(int));
                for (int ii = 0; ii < P; ii++){ 
                    bsp_put(ii, &empty, &empty, 0, sizeof(bool));
                }
            }
        }
        bsp_sync();

        if (empty == false){
            if (pid == diag_pid){
                if (winner_row_array[diag_pid] != -1){
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
            if (row_name == k % N && row_name != winner_row % N){
                tot_ones_row = find_ones_row(&head, all_ones_row, k);
                bsp_put(winner_row % N + (col_name % M) * N, all_ones_row, entries_received, 0, tot_ones_row * sizeof(int));
                delete_multiple_row(&head, k);
            }

            // Processors that handle the winner row send to those handling row k
            if (row_name == winner_row % N && row_name != k % N){
                tot_ones_row = find_ones_row(&head, all_ones_row, winner_row);
                bsp_put(k % N + (col_name % M) * N, all_ones_row, entries_received, 0, tot_ones_row * sizeof(int));
                //Delete all of your ones in that row
                delete_multiple_row(&head, winner_row);
            }
            bsp_sync();

            // Update the list of nodes
            if (row_name == k % N && row_name == winner_row % N){
                swap_two_rows(&head, k, winner_row);
            }
            else if (row_name == k % N){
                for (int i = 0; i < column_dim; i++){
                    if (entries_received[i] != -1){
                        append_to_list(&head, k, entries_received[i]);
                    }
                    else{
                        break;
                    }
                }
            }
            else if (row_name == winner_row % N){
                for (int i = 0; i < column_dim; i++){
                    if (entries_received[i] != -1){
                        append_to_list(&head, winner_row, entries_received[i]);
                    }
                    else{
                        break;
                    }
                }
            }

            // Reset vars
            for (int i = 0; i < column_dim; i++){
                all_ones_row[i] = -1;
            }

            // Row reduction steps! Get rid of ones in the same column as "col"
            // =========================================================================================================================================
            // Processors handling row "k" send to others in their processor column
            if (row_name == k % N){
                tot_ones_row = find_ones_row(&head, all_ones_row, k);
                for (int l = 0; l < N; l++){
                    bsp_put(l + (col_name % M) * N, all_ones_row, num_to_add, 0, tot_ones_row * sizeof(int));
                }
            }

            if (col_name == col % M){
                tot_ones_col = find_ones_col_except(&head, all_ones_col, col, k);
                for (int l = 0; l < M; l++){
                    bsp_put(row_name + (l % M) * N, all_ones_col, flag_array, 0, tot_ones_col * sizeof(int));
                }
            }
            bsp_sync();

            // Do the updates
            add_rows(&head, flag_array, num_to_add, row_dim, column_dim);

            // Reset vars
            for (int i = 0; i < P; i++){
                winner_row_array[i] = -1;
            }
            for (int i = 0; i < column_dim; i++){
                entries_received[i] = -1;
                num_to_add[i] = -1;
            }
            for (int i = 0; i < row_dim; i++){
                flag_array[i] = -1;
            }

            k++;

            bsp_sync();
        }
    }

    // Delete the registers used for communication between processors
    bsp_pop_reg(&empty);
    bsp_pop_reg(winner_row_array);
    bsp_pop_reg(&winner_row);
    bsp_pop_reg(flag_array);
    bsp_pop_reg(num_to_add);

    free(all_ones_row);
    free(all_ones_col);

    // All processors reconstruct the matrix, adding their 1's when they need to
    reconstruct_mat(&head);

    // Clean memory from all the child nodes
    struct node *child = head.next;
    head.next = NULL;
    while (child != NULL){
        struct node *new_child = child->next;
        free(child);
        child = new_child;
    }

    bsp_end();

    // Print the global matrix mat
    for(int i = 0; i < num_qubits; i++){
      for(int j = 0; j < 2 * num_qubits; j++){
        printf("%d,", mat[i][j]);
      }
      printf("\n");
    }

    // Write the output into a CSV file
    write_output();

    // Print the number of ones
    printf("%d\n", num_of_ones);
}

int main(int argc, char **argv){

    // Declare that "parallelalg1" is our parallel function
    bsp_init(sparsealg1, argc, argv);

    // Read the starting sparse stabilizer tableau 
    read_csv();

    // Read the starting stabilizer tableau and the vector of phases.
    sparsealg1();

    exit(EXIT_SUCCESS);
}
