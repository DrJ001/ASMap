#include <cctype>
#include <string>
#include <algorithm>
bool is_ok(char c) {
  // std::isdigit has undefined behaviour for argument values other than EOF
  // that are not representable as unsigned char. char is signed on the usual
  // platforms, so any byte >= 0x80 (for example in a non-ASCII marker name)
  // would be passed as a negative int. Convert before the call.
  return (!(std::isdigit(static_cast<unsigned char>(c)) || c=='.'));
}
bool is_number(std::string& s)
{
  return !s.empty() && std::find_if(s.begin(), s.end(), is_ok) == s.end();
}
