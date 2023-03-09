/*
* Copyright (C) 2013- Acho Arnold
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*
* Modifed March 2015 by Brett McClintock for use in the R package 'multimark'
*
*/

#include "matrix.h"
#include <R.h>

void matrix_print (Matrix *matrix)
{
  int i,j;
  Rprintf("\n");
  for(i = 0; i < matrix->row_size;i++)
  {
    Rprintf("\t\t");
    for(j = 0; j < matrix->col_size;j++)
    {
      Rprintf("%9.2f", matrix->matrix_entry[i][j]);
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}

void matrix_print_part(Matrix *matrix, int start_index)
{
  int j,i;
  for (i = 0; i < matrix->row_size; ++i)\
  {
    for (j = start_index; j < matrix->col_size; ++j)
    {
      Rprintf("\t\t%9.2f", matrix->matrix_entry[i][j]);
    }
    Rprintf("\n");
  }
}

/*Function to create an identity matrix */
Matrix * matrix_callalloc(int matrix_size)
{
  Matrix *result = matrix_alloc(matrix_size, matrix_size);
  int i,j;
  for (i = 0; i < matrix_size; i += 1)
  {
    for (j = 0; j < matrix_size; j += 1)
    {
      if (j == i)
      {
        result->matrix_entry[i][j] = 1;
      }
      else
      {
        result->matrix_entry[i][j] = 0;
      }
    }
  }
  return result;
}

Matrix *matrix_alloc(int row_size, int col_size)
{
  int j;
  Matrix *new_matrix = malloc(sizeof(Matrix));
  //Allocating memory for the new matrix structure
  new_matrix->row_size = row_size;
  new_matrix->col_size = col_size;
  new_matrix->matrix_entry = malloc( new_matrix->row_size *sizeof(double *));
  for(j = 0 ; j < new_matrix->row_size ; j++)
  {
    new_matrix->matrix_entry[j] = malloc( new_matrix->col_size*sizeof(double) );
  }
  return new_matrix;
}

/*Copies Matrix1 into matrix2 */
void matrix_copy(Matrix *matrix1, Matrix *matrix2)
{
  int i, j;
  for (i = 0; i < matrix1->row_size; i += 1)
  {
    for (j = 0; j < matrix1->col_size; j += 1)
    {
      matrix2->matrix_entry[i][j] = matrix1->matrix_entry[i][j];
    }
  }
}

void matrix_free( Matrix *matrix)
{
  int j;
  for(j = 0 ; j < matrix->row_size ; j++)
  {
    free(matrix->matrix_entry[j]);
  }
  free(matrix->matrix_entry);
  free(matrix);
}

/*Function which divides all row entries by the value of a the diagonal */
void row_divide(Matrix *matrix, int pivot)
{
  int j;
  double divisor = matrix->matrix_entry[pivot][pivot],
  result;
  for(j = pivot; j < matrix->col_size; j++)
  {
    result = (matrix->matrix_entry[pivot][j] / divisor);
    matrix->matrix_entry[pivot][j] = result;
  }
}

/*Function to carry out row operations*/
void row_operation(Matrix *multiplier_matrix,Matrix *matrix, int pivot, int row_index)
{
  int j;
  double multiplier = (matrix->matrix_entry[row_index][pivot] / matrix->matrix_entry[pivot][pivot]);
  //Loop which checks if matrix is provided to store the multiplier
  if(multiplier_matrix != NULL)
  {
    multiplier_matrix ->matrix_entry[row_index][pivot] = multiplier;
  }
  for(j=0; j < matrix->col_size; j++)
  {
    matrix->matrix_entry[row_index][j] -= multiplier * matrix->matrix_entry[pivot][j];
  }
}

void matrix_row_reduce( Matrix *matrix, int zero_control )
{
  int pivot, row_index;
  //double multiplier;
  for( pivot = 0; pivot < matrix->row_size ; pivot++)
  {
    error_zeros(matrix, zero_control); //Function checks if there are too many zeros in a single row
    if(  (matrix->matrix_entry[pivot][pivot] != 1) && (matrix->matrix_entry[pivot][pivot] != 0)  )
    {
      row_divide(matrix, pivot);
    }
    for (row_index = pivot+1; row_index < matrix->row_size; row_index++)
    {
      if (matrix->matrix_entry[pivot][pivot] != 0)
      {
        row_operation(NULL,matrix, pivot, row_index);
      }
    }
    for(row_index = pivot-1; row_index >=0; row_index --)
    {
      if (matrix->matrix_entry[pivot][pivot] != 0)
      {
      row_operation(NULL,matrix, pivot, row_index);
      }
    }
  }
}

void LU_decompose(Matrix *upper_triangular, Matrix *lower_triangular)
{
  int pivot, row_index;
  //double multiplier;
  for( pivot = 0; pivot < upper_triangular->row_size ; pivot++)
  {
    error_zeros(upper_triangular, upper_triangular->col_size); //Function checks if there are too many zeros in a single row
    for (row_index = pivot+1; row_index < upper_triangular->row_size; row_index++)
    {
      if ( upper_triangular->matrix_entry[pivot][pivot] != 0)
      {
        row_operation(lower_triangular,upper_triangular, pivot, row_index);
      }
    }
  }
}

void matrix_subtract(Matrix *result, Matrix *matrix1, Matrix *matrix2)
{
  int i, j;
  if ( !(matrix_equal_size(matrix1, matrix2 )) || \
  !(matrix_equal_size(matrix2, result)))
  {
    Rprintf("ERROR: The matrices you are trying to subtract have different sizes: please report to <brett.mcclintock@noaa.gov> \n");
    return;
  }
  for(i = 0; i < matrix1->row_size; i += 1)
  {
    for (j = 0; j < matrix1->col_size; j += 1)
    {
      result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] - matrix2->matrix_entry[i][j];
    }
  }
}

void matrix_add(Matrix *result, Matrix *matrix1, Matrix *matrix2)
{
  int i, j;
  if ( !(matrix_equal_size(matrix1, matrix2 )) || \
  !(matrix_equal_size(matrix2, result)))
  {
    Rprintf("ERROR: The matrices you are trying to add have different sizes: please report to <brett.mcclintock@noaa.gov> \n");
    return;
  }
  for (i = 0; i < matrix1->row_size; i += 1)
  {
    for (j = 0; j < matrix1->col_size; j += 1)
    {
    result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] + matrix2->matrix_entry[i][j];
    }
  }
}

void matrix_invert(Matrix *inverse_matrix)
{
  int j,k;
  /*Temporal matrix used in this function */
  Matrix *temp_matrix = matrix_alloc(inverse_matrix->row_size, inverse_matrix->col_size *2);
  matrix_copy(inverse_matrix, temp_matrix);
  /* Adding an identity matrix at the end of the temporal matrix */
  for(j = 0; j< temp_matrix->row_size; j++)
  {
    for(k = inverse_matrix->col_size; k < temp_matrix->col_size; k++)
    {
      if( (j+inverse_matrix->col_size) == k)
      {
        temp_matrix->matrix_entry[j][k] = 1;
      }
      else
      {
        temp_matrix->matrix_entry[j][k] = 0;
      }
    }
  }
  matrix_row_reduce(temp_matrix, temp_matrix->row_size);
  /*  for(j = 0; j< inverse_matrix->row_size; j++)
  {
    for(k = 0; k < inverse_matrix->col_size; k++)
    {
      Rprintf("matrix %d %d %f \n",j+1,k+1,inverse_matrix->matrix_entry[j][k]);
    }
    for(k = 0; k < temp_matrix->col_size; k++)
    {
      Rprintf("temp_matrix %d %d %f \n",j+1,k+1,temp_matrix->matrix_entry[j][k]);
    }
  }*/
  /* Copying the inverse matrix from the temp_matrix to the invse_matrix */
  for(j = 0; j< temp_matrix->row_size; j++)
  {
    for(k = inverse_matrix->col_size; k < temp_matrix->col_size; k++)
    {
      inverse_matrix->matrix_entry[j][k-inverse_matrix->col_size] = temp_matrix->matrix_entry[j][k];
      //Rprintf("matrix %d %d %f \n",j+1,k-inverse_matrix->col_size+1,temp_matrix->matrix_entry[j][k]);
    }
  }/*
      for(j = 0; j< inverse_matrix->row_size; j++)
  {
    for(k = 0; k < inverse_matrix->col_size; k++)
    {
      Rprintf("matrix %d %d %f \n",j+1,k+1,inverse_matrix->matrix_entry[j][k]);
    }
  }*/
  matrix_free(temp_matrix);
}

int matrix_equal_size( Matrix *matrix1, Matrix *matrix2)
{
  return (matrix1->row_size == matrix2->row_size && \
  matrix1->col_size == matrix2->col_size);
}

/*
This function checks if there is a line containing too many zero's and it exits
if such a line is found
*/
void error_zeros( Matrix *matrix, int control_index)
{
  int i,j,count;
  for(i=0; i<matrix->row_size; i++)
  {
    count=0;
    for(j = 0; j < matrix->col_size; j++)
    {
      if( matrix->matrix_entry[i][j] == 0)
      {
        count++;
      }
      else
      {
        return;
      }
      if(count == control_index)
      {
        Rprintf("Process fail because row %d contains %d zeros: please report to <brett.mcclintock@noaa.gov> \n",i+1,control_index);
        matrix_print(matrix);
        return;
      }
    }
  }
}
