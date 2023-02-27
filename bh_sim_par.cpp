#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include <mpi.h> 
#include "lib.h"

Vec2* position;
Vec2* i_velocity;
Vec2* velocity;
Vec2* force;
Cell* root_cell;
double* mass;
double* radius;


/* MPI datatypes for exchanging particles */
MPI_Datatype MPI_VEC2;


int N;               // particle count
int TIME;            
int rank;            
int size;            
int part_size;       // Number of particles each processor is responsible for
int pindex;          // The pindex points to the slot in the vectors/arrays that contains 
                     // data concerning the current processor, i.e pindex = (rank * part_size)


int name_length;                    // Length of processor name
char name[MPI_MAX_PROCESSOR_NAME];  // Buffer to hold processor name
double start_time;

double generate_rand(){
   return rand()/((double)RAND_MAX + 1);
}

void initialize(){


    for (int i = 0; i < N; i++) {
      mass[i] = M_UNKNOWN * generate_rand();
      radius[i] = RBOUND * generate_rand();
      position[i].x = generate_rand() * (X_WIDTH - RBOUND);
      position[i].y = generate_rand() * (Y_HEIGHT - RBOUND);
      i_velocity[i].x = 2 * generate_rand() - 1;
      i_velocity[i].y = 2 * generate_rand() - 1; 
   }
}


int check_collision(int index1, int index2) {

   if (pow((position[index1].x - position[index2].x), 2.0) + 
       pow((position[index1].y - position[index2].y), 2.0) <
       pow((radius[index1] + radius[index2]), 2.0)) {
       
       // Collision detected
       return 1;

   }   
   return 0;
}


double compute_distance(Vec2 a, Vec2 b){
    return sqrt(pow((a.x - b.x), 2.0) +
               pow((a.y - b.y), 2.0));
}

/*
 * Since the initial positions and radii were generated randomly
 * the system would have inherited a set of particles that were
 * already in collision. This method reinitializes the radius of
 * the particles such that the particles are not in collision with
 * each other.
 */
void reinitialize_radius() {

   int i, j;
   
   for (i = 0; i < N; i++) {

      for (j = i + 1; j < N; j++) {

         if (check_collision(i, j)) {
            double d = compute_distance(position[i], position[j]);
            radius[i] = radius[j] = d/2.0;
         }
      }
   }
}

void compute_force(){

   for (int i = 0; i < part_size; i++) {

      force[i] = Vec2();

      for (int j = 0; j < N; j++){

         if (j == (i + pindex)) continue; // avoid computation for same bodies

         double d = compute_distance(position[i + pindex], position[j]);
         double f = (G * (mass[i + pindex] * mass[j]) / (pow(d, 2.0)));

         force[i] = force[i] + ((position[j] - position[i + pindex]) / d) * f;
      }
   }
}

void compute_velocity(){
   for (int i = 0; i < part_size; i++) {
      velocity[i] = velocity[i] + (force[i] / mass[i + pindex]) * DELTAT;
   }
}

void compute_positions(){

   for (int i = 0; i < part_size; i++) {
      position[i + pindex] += velocity[i] * DELTAT;

      // Check if particles attempt to cross boundary      
      if ((position[i + pindex].x + radius[i + pindex]) >= X_WIDTH ||
          (position[i + pindex].x - radius[i + pindex]) <= 0)
         velocity[i].x *= -1;
      else if ((position[i + pindex].y + radius[i + pindex] >= Y_HEIGHT) || 
               (position[i + pindex].y - radius[i + pindex]) <= 0)
         velocity[i].y *= -1;     
   }
}

Cell* create_cell(double x, double y){
    Cell* cell = (Cell*)malloc(sizeof(Cell));
    cell->mass = 0;
    cell->no_childs = 0;
    cell->index = -1;
    cell->cx = 0;
    cell->cy = 0;
    cell->width = x;
    cell->height = y;  
    return cell;
}

void set_children(Cell* cell, double width, double heigth){

   // Set location of new cells
   cell->children[0]->x = cell->x;
   cell->children[0]->y = cell->y;

   cell->children[1]->x = cell->x + width;
   cell->children[1]->y = cell->y;

   cell->children[2]->x = cell->x;
   cell->children[2]->y = cell->y + heigth;

   cell->children[3]->x = cell->x + width;
   cell->children[3]->y = cell->y + heigth;

}

void generate_childs(Cell* cell) {
   
   double width  = cell->width / 2.0;
   double height = cell->height / 2.0;

   cell->no_childs = 4;   
   
   // Create and initialize new subcells   
   for (int i = 0; i < cell->no_childs; i++) {
      cell->children[i] = create_cell(width, height);
   }   
   set_children(cell, width, height);   
}

int locate_child(Cell* cell, int index){

    if (position[index].x > cell->children[3]->x){
        if (position[index].y > cell->children[3]->y){
            return 3;
        }
        else{
            return 1;
        }
    }
    else{
        if (position[index].y > cell->children[3]->y){
            return 2;
        }
        else{
            return 0;
        }      
    }
}


void add_to_cell(Cell* cell, int index) {

   if (cell->index == -1) {         
      cell->index = index;
      return;         
   }         
   generate_childs(cell);

   int sc1 = locate_child(cell, cell->index);
   cell->children[sc1]->index = cell->index;   

   int sc2 = locate_child(cell, index);

   if (sc1 == sc2)
      add_to_cell(cell->children[sc1], index);
   else 
      cell->children[sc2]->index = index;  
}


void tree_generator() {
   
   // Initialize root of Quadtree
   root_cell = create_cell(X_WIDTH, Y_HEIGHT);
   root_cell->index = 0;
   root_cell->x = 0;
   root_cell->y = 0;

   for (int i = 1; i < N; i++) {

      Cell* cell = root_cell;
      while (cell->no_childs != 0){
         int sc = locate_child(cell, i);
         cell = cell->children[sc];
      }      
      add_to_cell(cell, i);
   }
}


Cell* compute_cell_prop(Cell* cell){
   
   if (cell->no_childs == 0) {
      if (cell->index != -1){
         cell->mass = mass[cell->index];
         return cell;
      }
   }
   else {      
      Vec2 t = Vec2();
      for (int i = 0; i < cell->no_childs; i++) {
         Cell* temp = compute_cell_prop(cell->children[i]);
         if (temp != NULL) {
            cell->mass += temp->mass;
            t += position[temp->index] * temp->mass;           
         }
      }      
      // Compute center of mass
      cell->cx = t.x / cell->mass;
      cell->cy = t.y / cell->mass;
   
      return cell;
   }
   return NULL;
}

void compute_force_cell(Cell* cell, int index) {
   double d = compute_distance(position[index], position[cell->index]);

   // Compute grativational force according to Newtonian's law
   double f = (G * (mass[index] * mass[cell->index]) / (pow(d, 2.0)));

   // Resolve forces in each direction
   force[index - pindex] += ((position[cell->index] - position[index]) / d) * f;   
}

void compute_force_tree(Cell* cell, int index) {
   
   if (cell->no_childs == 0) {
      if (cell->index != -1 && cell->index != index) {
         compute_force_cell(cell, index);
      }
   }
   else {
      double d = compute_distance(position[index], position[cell->index]);
      
      if (THETA > (cell->width / d)){ 
         // Use approximation
         compute_force_cell(cell, index);         
      }
      else {
         int i;
         for (i = 0; i < cell->no_childs; i++) {
            compute_force_tree(cell->children[i], index);
         }
      }      
   }
}

void compute_force_space(){

   int i, j;   

   for (i = 0; i < part_size; i++) {

      force[i] = Vec2();
      compute_force_tree(root_cell, i + pindex);
   }
}


void delete_tree(Cell* cell) {
   
   if (cell->no_childs == 0) {
      free(cell);
      return;
   }

   int i;
   for (i = 0; i < cell->no_childs; i++) {
      delete_tree(cell->children[i]);
   }

   free(cell);
}

void init_velocity(){
   int i;
   for (i = 0; i < part_size; i++){
      velocity[i] = Vec2();
   }
}

void write_positions() {

   FILE* file;

   file = fopen("pdist.data", "w");

   if (file == NULL) {
       fprintf(stderr,"Cannot open output file\n");
       exit (0);
   }

   int i;
   for (i = 0; i < N; i++) {
       fprintf(file, "%f, %f\n", position[i].x, position[i].y);
   }

   fclose(file);
}


void run_simulation(){

   // Broadcast mass and position to all members in the group
   MPI_Bcast(mass, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(position, N, MPI_VEC2, 0, MPI_COMM_WORLD);
   MPI_Scatter(i_velocity, part_size, MPI_VEC2, velocity, part_size, MPI_VEC2, 0, MPI_COMM_WORLD);

   for (int i = 0; i < TIME; i++) {
      tree_generator();
      compute_cell_prop(root_cell);
      compute_force_space();
      delete_tree(root_cell);

      compute_velocity();
      compute_positions();
      MPI_Allgather(position + (rank * part_size), part_size, MPI_VEC2, 
                    position, part_size, MPI_VEC2, MPI_COMM_WORLD); 
   }

   if (rank == 0){
      write_positions();
      double end_time =  MPI_Wtime();
      std::cout << "N = " << N << ", IT ="<< TIME;
      std::cout << ", Time elapsed = " << end_time - start_time<< " seconds.\n";

      FILE* file;
      file = fopen("speed_par.txt", "a");
      fprintf(file, "%d, %d, %f\n", N, TIME, end_time - start_time);

      fclose(file);
   }
}

int main(int argc, char* argv[]){

   MPI_Init(&argc, &argv);  
   // Initialise problem parameters
   if (argc >= 2) 
      sscanf(argv[1], "%i%", &N);
   else
      N = DEFAULT_N;

   if (argc >= 3)
      sscanf(argv[2], "%i%", &TIME);   
   else
      TIME = DEFAULT_TIME;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   // Create and commit new MPI Types
   MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_VEC2);
   MPI_Type_commit(&MPI_VEC2);

   // Identify processor
   MPI_Get_processor_name(name, &name_length);
   // printf("Rank=%d, Processor=%s, N=%d, TIME=%d\n", rank, name, N, TIME);
   start_time = MPI_Wtime();

   part_size = N / size;
   pindex = rank * part_size;

   mass = (double *) malloc(N * sizeof(double));
   radius = (double *) malloc(N * sizeof(double));
   position = (Vec2 *) malloc(N * sizeof(Vec2));
   i_velocity = (Vec2 *) malloc(N * sizeof(Vec2));
   velocity = (Vec2 *) malloc(part_size * sizeof(Vec2));
   force = (Vec2 *) malloc(part_size * sizeof(Vec2));

   // Initialize velocity array for each process
   init_velocity();

   // Let the master initialize the space
   if (rank == 0){
      initialize();
   }
      
   // Run the N-body simulation
   run_simulation();

   // Terminate MPI execution env.
   MPI_Finalize();

   return 0;                   
}