/* M.Nielsen July 2008, mniel@cbs.dtu.dk */

#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>

extern int add_one(int i); //declare file



main(int argc, char *argv[])
{

      //printf( "Hello World\n" );
      //int i;
      
      //for(i=0;i<=100;i++){
      //      printf("%d \t",i);
      //      printf("%d \n",add_one(i));
      //}


      FILE *ifp; //declare variable for reading input file.
      char *fname = "1A68_HUMAN.sprot";
      char line[1024];
      char dummy[256], id[256], seq[1024];
      int i=0;
      int j=0;
      

      ifp = fopen(fname, "r"); //open file


      if(ifp == NULL)
      {
            printf("can't open file\n");
            return;
      }      

      while( fgets(line,sizeof line,ifp) != NULL) 
      {
            sscanf( line, "%s %s", dummy, id);

            if(i==0)
            {                  
                  if(strcmp(dummy,"ID")==0)
                  {
                        printf(">");
                        printf("%s\n", id);
                  }
                  else if(strcmp(dummy,"SQ")==0)    
                        i=1;
            }
            else
            {
                  if(strcmp(dummy,"//") != 0 )
                  {
                        for(j=0;j<strlen(line);j++)
                        {
                              char t[64];
                              strncpy(t, line + j, 1 );
                              if(isspace(line[j])==0)
                                    strcat(seq,t);
                        }
                        strcat(seq,"\n");
                  }
                  else
                  {
                        printf("%s",seq);           
                  }
                  
            }
      }
      fclose( ifp );
}


      
      

int add_one( int i )
{
      int j = i+1;
      return(j);     
}
