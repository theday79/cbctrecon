/****************************************************************************
** Modified from qcombobox.h from the Qt Project by Andreas Gravgaard Andersen
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtWidgets module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QCOMBOBOX_H
#define QCOMBOBOX_H
#include <qglobal.h>
#include <qicon.h>
#include <qvalidator.h>

#include <QRegularExpression>

class QAbstractItemView;
class QLineEdit;
class QCompleter;

/** determine if the pattern matches the string using Qt::MatchFlags
    \param str the string
    \param pattern the pattern to find
    \param flags any combination of the follow Qt flags
                - Qt::MatchFixedString
                - Qt::MatchContains
                - Qt::MatchStartsWith
                - Qt::MatchEndsWith
                - Qt::MatchRegExp (overrides all flags above)
                - Qt::MatchCaseSensitive
    \returns true if the pattern is found in the string
*/
static bool
QString_Matches(const QString &str, const QString &pattern,
                const Qt::MatchFlags &flags = (Qt::MatchCaseSensitive |
                                               Qt::MatchFixedString)) {
  if (flags.testFlag(Qt::MatchRegExp)) {
    QRegularExpression::PatternOptions options =
        QRegularExpression::NoPatternOption;
    if (!flags.testFlag(Qt::MatchCaseSensitive)) {
      options = QRegularExpression::CaseInsensitiveOption;
    }
    QRegularExpression regex(pattern, options);
    return regex.match(str).hasMatch();
  }

  auto cs = Qt::CaseSensitive;
  if (!flags.testFlag(Qt::MatchCaseSensitive)) {
    cs = Qt::CaseInsensitive;
  }
  if (flags.testFlag(Qt::MatchContains)) {
    return str.contains(pattern, cs);
  }
  if (flags.testFlag(Qt::MatchStartsWith)) {
    if (str.startsWith(pattern, cs)) {
      return true;
    }
  }
  if (flags.testFlag(Qt::MatchEndsWith)) {
    if (str.endsWith(pattern, cs)) {
      return true;
    }
  }
  if (flags.testFlag(Qt::MatchFixedString)) {
    return str.compare(pattern, cs) == 0;
  }
  return false;
}

class MyQComboBox {

public:
  MyQComboBox() = default;
  ~MyQComboBox() = default;

  int count() const { return m_container.size(); }
  void setMaxCount(const int max) { m_container.reserve(max); }
  int maxCount() const { return m_container.capacity(); }

  int findText(const QString &text,
               const Qt::MatchFlags flags = static_cast<Qt::MatchFlags>(
                   Qt::MatchExactly | Qt::MatchCaseSensitive)) const {
    return findData(text, Qt::DisplayRole, flags);
  }
  int findData(const QVariant &data, int role = Qt::UserRole,
               const Qt::MatchFlags flags = static_cast<Qt::MatchFlags>(
                   Qt::MatchExactly | Qt::MatchCaseSensitive)) const {
    auto i_out = 0;
    for (auto &item : m_container) {
      const auto match =
          QString_Matches(item.first, data.toString(), flags) ||
          QString_Matches(item.second.toString(), data.toString(), flags);
      if (match) {
        return i_out;
      }
      ++i_out;
    }
    return -1;
  }

  enum InsertPolicy {
    NoInsert,
    InsertAtTop,
    InsertAtCurrent,
    InsertAtBottom,
    InsertAfterCurrent,
    InsertBeforeCurrent,
    InsertAlphabetically
  };

  InsertPolicy insertPolicy() const { return m_insert_policy; }
  void setInsertPolicy(const InsertPolicy policy) { m_insert_policy = policy; }

  bool isEditable() const { return m_editable; }
  void setEditable(const bool editable) { m_editable = editable; }

  int currentIndex() const { return m_current_index; }
  QString currentText() const { return m_container.at(m_current_index).first; }

  QVariant currentData(int role = Qt::UserRole) const {
    return m_container.at(m_current_index).second;
  }

  QString itemText(const int index) const {
    return m_container.at(index).first;
  }

  QVariant itemData(const int index, int role = Qt::UserRole) const {
    return m_container.at(index).second;
  }

  void addItem(const QString &text, const QVariant &userData = QVariant()) {
    insertItem(count(), text, userData);
  }

  void addItem(const QIcon &icon, const QString &text,
               const QVariant &userData = QVariant()) {
    insertItem(count(), icon, text, userData);
  }

  void addItems(const QStringList &texts) { insertItems(count(), texts); }

  void insertItem(const int index, const QString &text,
                  const QVariant &userData = QVariant()) {
    insertItem(index, QIcon(), text, userData);
  }

  void insertItem(const int index, const QIcon &icon, const QString &text,
                  const QVariant &userData = QVariant()) {
    if (!m_editable) {
      return;
    }
    QVariant data;
    if (!userData.isValid()) {
      data = QVariant(text);
    } else {
      data = userData;
    }
    const auto pair = std::make_pair(text, data);
    const auto it = m_container.begin();
    m_container.insert(it + index, pair);
  }

  void insertItems(const int index, const QStringList &texts) {
    size_t i = index;
    for (auto &item : texts) {
      insertItem(i, item);
      ++i;
    }
  }

  void removeItem(const int index) {
    const auto it = m_container.begin();
    m_container.erase(it + index);
  }

  void setItemText(const int index, const QString &text) {
    if (!m_editable || index < 0 ||
        m_container.size() <= static_cast<size_t>(index)) {
      return;
    }
    m_container.at(index).first = text;
  }

  void setItemData(const int index, const QVariant &value,
                   int role = Qt::UserRole) {
    if (!m_editable || index < 0 ||
        m_container.size() <= static_cast<size_t>(index)) {
      return;
    }
    m_container.at(index).second = value;
  }

  void clear() {
    m_container.clear();
    m_current_index = 0;
  }

  void setCurrentIndex(const int index) { m_current_index = index; }

  void setCurrentText(const QString &text) {
    m_container.at(m_current_index) = std::make_pair(text, QVariant(text));
  }

private:
  Q_DISABLE_COPY(MyQComboBox)
  std::vector<std::pair<QString, QVariant>> m_container;
  size_t m_current_index = 0;
  bool m_editable = true;
  InsertPolicy m_insert_policy = InsertAtBottom;
};

QT_END_NAMESPACE

#endif // QCOMBOBOX_H
