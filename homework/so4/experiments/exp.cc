#include "exp.hpp"

#include <fstream>
#include <sstream>

Input *DataGen::readFile(string filename)
{
    fstream file(filename, ios::in);
    Input *input = new Input();

    if (file.is_open())
    {
        string line;
        int i = -1;

        while (getline(file, line))
        {
            istringstream iss(line);

            if (i == -1)
            {
                iss >> input->n;

                input->data = new data_t[input->n];
            }
            else
            {
                istringstream iss(line);

                iss >> input->data[i];
            }
            ++i;
        }

        file.close();
    }

    return input;
}

void DataGen::genRand(int n, data_t *&data)
{
    for (int i = 0; i < n; ++i)
    {
        data[i] = static_cast<data_t>(i);
    }
}

void DataGen::genSequential(int n, data_t *&data)
{
    for (int i = 0; i < n; ++i)
    {
        data[i] = static_cast<data_t>(rand() % 10000);
    }
}
