#include <utilities/strings.hpp>

#include <cstring>
#include <stdio.h>

#include <iostream>
#include <algorithm>

// String helpers
bool endsWith(std::string_view str, std::string_view suffix) {
    return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

bool startsWith(std::string_view str, std::string_view prefix) {
    return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}

void strip(std::string& str, char c) { 
    str.erase(remove(str.begin(), str.end(), c), str.end()); 
}

std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token     = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

int substrSscanf(const char* src, int start, int end, const char* format, void* result) {
    char* substr = (char*)malloc(end - start + 2);
    memcpy(substr, &src[start], end - start + 1);
    substr[end - start + 1] = '\0';

    int matches = sscanf(substr, format, result);

    free(substr);
    return matches;
}

// @note relies on end - start + 2 minimum chars allocated for result
void substrString(const char* src, int start, int end, char* result) {
    memcpy(result, &src[start], end - start + 1);
    result[end - start + 1] = '\0';
}

// @note relies on end - start + 2 minimum chars allocated for result
void substrChar(const char* src, int pos, char* result) { (*result) = src[pos]; }

void substrInt(const char* src, int start, int end, int* result) { substrSscanf(src, start, end, " %d", result); }

void substrFloat(const char* src, int start, int end, float* result) { substrSscanf(src, start, end, " %f", result); }