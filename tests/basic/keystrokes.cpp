#include <stdio.h>
#include <stdlib.h>  // for atexit()
#include <time.h>    // For measuring time

#if defined(__linux__) or defined(__unix__)
#include <fcntl.h>   // file control
#include <pthread.h> // Multithreading
#include <termios.h> // For changing terminal mode
#include <unistd.h>  // For changing terminal mode
#endif 

#include "../../basic.hpp"


#if defined(__linux__) or defined(__unix__)

static termios original;   // A struct to save the original state of terminal
static int process_reading_is_finished = 0; // For thread communication

void restore_terminal_mode();
void change_terminal_mode();
void *process_read_keyboard( void* );
void *process_print( void* );

/// Reads keyboard input
void *process_read_keyboard( void* ) {
  
  change_terminal_mode(); 
  
  const auto t_start = clock();
  
  while( clock() - t_start < 10 * CLOCKS_PER_SEC ) {
    
    char c;
    int bytesRead = read( STDIN_FILENO, &c, 1 );
    
    if( bytesRead > 0 ) {

      if( c == 27 ) {
        printf("You pressed ESC\n");
        break;
      } else {
        printf("You pressed %c\n", c );
      }

    }

  }
  
  process_reading_is_finished = 1;
  
  return nullptr;
}

/// Doing Stuff while listening to keyboard
void *process_print( void* ) {
  
  const auto t_start = clock();
  
  while( !process_reading_is_finished && ( clock() - t_start ) < 10 * CLOCKS_PER_SEC ) { 
    sleep(1);
    LOG << "I am Printing! Time: " << ( clock() - t_start ) << nl;
  }
  
  LOG << "Printing Thread Finished!\n";
  
  display_mallinfo(false);
  
  return nullptr;

}

/// This function changes the terminal mode.
void change_terminal_mode() {
  struct termios raw;
  
  tcgetattr(STDIN_FILENO, &raw); 
  // Save the state of the terminal to struct raw
  // STDIN_FILENO is from <stdlib.h>
  // tcgetattr() from <termios.h>
  
  tcgetattr(STDIN_FILENO, &original);
  
  atexit(&restore_terminal_mode);
  // Revert to canonical mode when exiting the program
  // atext() from <stdlib.h>
  
  raw.c_lflag &= ~(ECHO | ICANON); 
  raw.c_cc[VMIN]  = 1;
  raw.c_cc[VTIME] = 0;
  // Turn off canonical mode
  // Turn off ECHO mode so that keyboard is not
  // printing to terminal
  // ICANON and ECHO is bitflag. ~ is binary NOT operator

  tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
  // Set the terminal to be in raw mode
  // tcsetattr() from <termios.h>

  int fcntlFlags = fcntl( STDIN_FILENO, F_GETFL, 0 );
  fcntl( STDIN_FILENO, F_SETFL, fcntlFlags | O_NONBLOCK );
}

void restore_terminal_mode() {
  tcsetattr(STDIN_FILENO, TCSAFLUSH, &original); 
  // Set terminal to original state
}

#endif 



int main() {

#if defined(__linux__) or defined(__unix__)

  // Start Multithreading
  pthread_t id_print, id_read;

  pthread_create(&id_print, nullptr, process_print, nullptr);
  pthread_create(&id_read, nullptr, process_read_keyboard, nullptr);

  pthread_join(id_print, nullptr);
  pthread_join(id_read, nullptr);

  
#endif 

  return 0;
}



