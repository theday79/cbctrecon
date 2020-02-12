#ifndef QTWRAP_H
#define QTWRAP_H

#include <filesystem>
#include <string>
#include <vector>

#include <QString>
#include <QstringList>

namespace fs = std::filesystem;

template <typename T, typename std::enable_if_t<
                          std::is_constructible_v<QString, T>, int> = 0>
QString to_qstr(const T &path) {
  return QString(path);
}

template <typename T,
          typename std::enable_if_t<std::is_same_v<T, fs::path>, int> = 0>
[[nodiscard]] QString to_qstr(const T &path) {
  return QString(path.string().c_str());
}

template <typename T,
          typename std::enable_if_t<std::is_same_v<T, std::string>, int> = 0>
[[nodiscard]] QString to_qstr(const T &path) {
  return QString(path.c_str());
}

template <typename T, typename std::enable_if_t<
                          std::is_constructible_v<fs::path, T>, int> = 0>
[[nodiscard]] fs::path to_path(const T &path) {
  return fs::path(path);
}

template <typename T,
          typename std::enable_if_t<std::is_same_v<T, QString>, int> = 0>
[[nodiscard]] fs::path to_path(const T &path) {
  return fs::path(path.toStdString());
}

template <typename T>
[[nodiscard]] std::vector<T> qstr_list_to_t_vec(const QStringList &qstr_list) {
  auto t_vec = std::vector<T>();
  std::transform(qstr_list.cbegin(), qstr_list.cend(),
                 std::back_inserter(t_vec),
                 [](const QString &qstr) { return T(qstr.toStdString()); });
  return t_vec;
}

#endif // QTWRAP_H
