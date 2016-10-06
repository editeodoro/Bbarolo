/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 BBarolo is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#ifndef ConsoleStream_H
#define ConsoleStream_H

// -- stl stuff
#include <iostream>
#include <streambuf>
#include <string>

// -- QT stuff
#include <QtGui>

class ConsoleStream : public std::basic_streambuf<char> {
    public:
        ConsoleStream(std::ostream *stream, QTextEdit* textEdit) {
            previousBuffer = NULL;
            logTextEdit = NULL;
            myStream = NULL;
            init(stream, textEdit);
        }

        ConsoleStream() {
            previousBuffer = NULL;
            logTextEdit = NULL;
            myStream = NULL;
        }

        ~ConsoleStream() {
            if (!myString.empty())
                logTextEdit->append(myString.c_str());
            free();
        }

       void setStream(std::ostream *stream) {
            free();
            myStream = stream;
            previousBuffer = stream->rdbuf();
            stream->rdbuf(this);
        }

        void setTextEdit(QTextEdit* text_edit) {
            logTextEdit = text_edit;
        }

        void init(std::ostream *stream, QTextEdit* textEdit) {
            setTextEdit(textEdit);
            setStream(stream);
        }

        void free() {
            if (previousBuffer != NULL && myStream != NULL)
                myStream->rdbuf(previousBuffer);
        }

    protected:
        virtual int_type overflow(int_type v) {
            if (v == '\n') {
                logTextEdit->append(myString.c_str());
                myString.erase(myString.begin(), myString.end());
            } else myString += v;

            return v;
        }

        virtual std::streamsize xsputn(const char *p, std::streamsize n) {
            myString.append(p, p + n);

            std::string::size_type pos = 0;
            while (pos != std::string::npos) {
                pos = myString.find('\n');
                if (pos != std::string::npos) {
                    std::string tmp(myString.begin(), myString.begin() + pos);
                    logTextEdit->append(tmp.c_str());
                    myString.erase(myString.begin(), myString.begin() + pos + 1);
                }
            }

            return n;
        }

    private:
        std::ostream *myStream;
        std::streambuf *previousBuffer;
        std::string myString;
        QTextEdit* logTextEdit;
};

 #endif
