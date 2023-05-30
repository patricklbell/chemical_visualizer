#ifndef UTILITIES_STRING_HPP
#define UTILITIES_STRING_HPP

#include <string>
#include <vector>

// String helpers
bool startsWith(std::string_view str, std::string_view prefix);
bool endsWith(std::string_view str, std::string_view suffix);
void strip(std::string &str, char c);
std::vector<std::string> split(std::string s, std::string delimiter);

int substrSscanf(const char *src, int start, int end, const char *format, void *result);
void substrString(const char *src, int start, int end, char *result);
void substrChar(const char *src, int pos, char *result);
void substrInt(const char *src, int start, int end, int *result);
void substrFloat(const char *src, int start, int end, float *result);

#endif    // UTILITIES_STRING_HPP