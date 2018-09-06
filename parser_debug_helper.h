#ifndef PARSER_DEBUG_HELPER
#define PARSER_DEBUG_HELPER

#include <stdio.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
	# define __BEGIN_DECLS extern "C" {
	# define __END_DECLS }
#else
	# define __BEGIN_DECLS /* empty */
	# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void printf_red(char *str){
	printf(ANSI_COLOR_RED "%s\t\t", str);
	printf(ANSI_COLOR_RESET);
}
void printf_green(char *str){
	printf(ANSI_COLOR_GREEN "%s\t\t", str);
	printf(ANSI_COLOR_RESET);
}
void printf_blue(char *str){
	printf(ANSI_COLOR_BLUE "%s\t\t", str);
	printf(ANSI_COLOR_RESET);
}
void printf_yellow(char *str){
	printf(ANSI_COLOR_YELLOW "%s\t\t", str);
	printf(ANSI_COLOR_RESET);
}
void printf_magenta(char *str){
	printf(ANSI_COLOR_MAGENTA "%s\t\t", str);
	printf(ANSI_COLOR_RESET);
}
void printf_cyan(char *str){
	printf(ANSI_COLOR_CYAN "%s\t\t", str);
	printf(ANSI_COLOR_RESET);
}

__END_DECLS

#endif 