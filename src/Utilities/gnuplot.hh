//------------------------------------------------------------------
// gnuplot.hh: Some classes in order to interact with Gnuplot
//------------------------------------------------------------------

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
    Internet email: enrico.diteodoro@gmail.com
-----------------------------------------------------------------------*/

#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <cstdio>
#include <string>
#include <cstdarg>

class Gnuplot_Base 
{
protected:
  
    std::string exec;
    std::string cmdline;
    FILE * fp;
  
public:
    explicit Gnuplot_Base(const char * execname=NULL) : fp(NULL) {
        if(execname) {
            exec = execname;
        }
        else {
            exec = "gnuplot";
        }
    }
  
    virtual ~Gnuplot_Base() {
        end();
    }

    enum { GNUPLOT_OK = 0,
           GNUPLOT_NG,
           GNUPLOT_ALREADY_OPEN,
           GNUPLOT_NOT_OPEN,
           GNUPLOT_CANNOT_EXEC
    };
  
    virtual int begin(const char * cmd_opt=NULL, const char * reserved=NULL) = 0;
    virtual int end() {return 0;}

    int flush() {
        if(fp)::fflush(fp); 
        return GNUPLOT_OK;
    }

    int begin_data() {
        if(!fp) return GNUPLOT_NOT_OPEN;
        return flush();
    }

    int end_data() {
        if( ! fp ) return GNUPLOT_NOT_OPEN;
        ::fprintf(fp, "e\n");
        return flush();
    }
  
    int command(const char * cmd, ...) {
        if(!fp) return GNUPLOT_NOT_OPEN;
        if(!cmd) return GNUPLOT_OK;
        {va_list ap; va_start(ap, cmd); ::vfprintf(fp, cmd, ap); va_end(ap);}
        return GNUPLOT_OK;
    }

    int commandln(const char * cmd, ...) {
        if(!fp) return GNUPLOT_NOT_OPEN;
        if(!cmd) return GNUPLOT_OK; 
        {va_list ap; va_start(ap, cmd); ::vfprintf(fp, cmd, ap); va_end(ap);}
        ::fprintf(fp, "\n");
        return GNUPLOT_OK;
    }
};


class Gnuplot_Pipe : public Gnuplot_Base
{
public:
    virtual int begin(const char * cmd_opt=NULL, const char * reserved=NULL) {
        if(fp) return GNUPLOT_ALREADY_OPEN;
        cmdline = exec;
        if(cmd_opt) {
            cmdline += " ";
            cmdline += cmd_opt;
        }
        fp = ::popen( cmdline.c_str(), "w" );
        if(!fp) return GNUPLOT_CANNOT_EXEC;
        return GNUPLOT_OK;
    }
  
    virtual int end() {
        if(!fp) return GNUPLOT_OK; 
        ::fflush(fp);    
        ::pclose(fp);
        fp = NULL;
        return GNUPLOT_OK;
    }
  
};


class Gnuplot_Tmpfile : public Gnuplot_Base
{
protected:
    std::string tmpfilename;
    bool remove_file;
  
public:
    virtual int begin(const char * cmd_opt=NULL, const char * filename=NULL) {
        if(fp) return GNUPLOT_ALREADY_OPEN;
        cmdline = exec;
        if(cmd_opt) {
            cmdline += " ";
            cmdline += cmd_opt;
        }
        if(filename) {
            fp = ::fopen(filename, "w");
            tmpfilename = filename;
            remove_file = false;
        }
        else {
            char buf[] = "/tmp/gnuplot_XXXXXX";
            int fd = ::mkstemp(buf);
            if(fd==-1) return GNUPLOT_NG;
            fp = ::fdopen( fd, "w" );
            tmpfilename = buf;
            remove_file = true;
        }
    
        if(!fp) {
            tmpfilename.clear();
            return GNUPLOT_CANNOT_EXEC;
        }
        
        return GNUPLOT_OK;
    }
  
    virtual int end() {
        if(!fp) return GNUPLOT_OK;
        ::fflush(fp);
        ::fclose(fp);
        fp = NULL;
        if(! ::system((cmdline + " " + tmpfilename).c_str())) 
            return GNUPLOT_CANNOT_EXEC;
        if(remove_file) ::remove( tmpfilename.c_str());
        return GNUPLOT_OK;
    }
};

// use pipe version as default
typedef Gnuplot_Pipe Gnuplot;

#endif 
