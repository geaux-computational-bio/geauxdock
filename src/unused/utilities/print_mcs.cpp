#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <cstdlib>
#include <cstdio>

using namespace std;


void
LoadLHM (std::string fnpath)
{
    ifstream ifn (fnpath.c_str ());
    if (!ifn.is_open ()) {
        cout << "Failed to open " << fnpath << endl;
        exit (EXIT_FAILURE);
    }

    std::string line;
    while (getline (ifn, line)) {
        //cout << line.length () << endl;
        //if (line.substr (0, 3) == "KDE") ...

        std::istringstream iss (line);
        std::vector < std::string > tokens;
        while (iss) {
            std::string t;
            iss >> t;
            tokens.push_back (t);
        }

        if (tokens.size () <= 1)
            continue;

        if (tokens[0] == "MCS") {
            //float tcc = atof (tokens[2].c_str ());
            int num = atoi (tokens[3].c_str ());

            //cout << line << endl;

            printf ("valid mcs points %3d:    ", num);
            for (int i = 0; i < num; ++i) {
                const int offset = 4;
                int loc = atoi (tokens[4 * i + 0 + offset].c_str ());
                //float x = atof (tokens[4 * i + 1 + offset].c_str ());
                //float y = atof (tokens[4 * i + 2 + offset].c_str ());
                //float z = atof (tokens[4 * i + 3 + offset].c_str ());
                printf ("%2d ", loc);
            }
            printf ("\n");
        }

    }

    ifn.close ();
}


int
main (int argc, char **argv)
{
    if (argc == 1) {
        printf ("%s <LHMFILE(*.ff)>\n", argv[0]);
        exit (EXIT_FAILURE);
    }

    LoadLHM (argv[1]);

    return 0;
}

