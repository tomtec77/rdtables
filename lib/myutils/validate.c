/**
 * @file validate.c
 * @brief Prompt user to authorize continuation of the program
 * @author Tomas E. Tecce
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "myutils.h"

unsigned int validate_continuation(void)
{
  int valid_input;           /* When 1, data is valid and loop is exited */
  char  user_input;  /* Handles user input, single character menu choice */

  valid_input = 0;
  while( valid_input == 0 ) {
    printf("Continue anyway (Y/N)?\n");
    (void)scanf("  %c", &user_input );
    user_input = toupper( user_input );
    if((user_input == 'Y') || (user_input == 'N') )  valid_input = 1;
    else  printf("\007Error: Invalid choice\n");
  }

  if (user_input == 'Y') return 1;
  else return 0;
}
