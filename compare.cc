/*************************************************************************
 *
 * Copyright (C) 2012. Vegard Nossum (vegardno@ifi.uio.no)
 * Copyright (C) 2010. Mladen Nikolic (nikolic@matf.bg.ac.rs)
 *
 * This program uses ALGLIB (www.alglib.net)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (www.fsf.org); either version 2 of the 
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public License is available at
 * http://www.fsf.org/licensing/licenses
 *
 *************************************************************************/

#include <fstream>
#include <iostream>

#include "GehanWilcoxon.hpp"

Sample read_sample(const char *filename)
{
	Sample s;
	std::ifstream in(filename);

	while (true) {
		std::string token;
		in >> token;
		if (!in.good())
			break;

		bool known = true;
		if (token == ">") {
			known = false;

			in >> token;
			if (!in.good())
				break;
		}

		s.addDatum(Datum(std::stod(token), known));
	}

	return s;
}

int main(int argc, char *argv[])
{
	if (argc != 3)
		exit(1);

	Sample s1 = read_sample(argv[1]);
	Sample s2 = read_sample(argv[2]);

	GehanWilcoxon test(s1, s2, 0);
	double rv = test.rValue();
	double rv2 = rv * rv;
	double zsum = 0.5 * log((1 + rv) / (1 - rv));
	double zvar = test.variance() / ((1 - rv2) * (1 - rv2));

	std::cout << "n1 = " << s1.size() << std::endl;
	std::cout << "n2 = " << s2.size() << std::endl;
	std::cout << "z = " << zsum / sqrt(zvar) << std::endl;
	std::cout << "r = " << rv << std::endl;
	std::cout << "p = " << 2 - 2 * normaldistribution(fabs(zsum / sqrt(zvar))) << std::endl;
	std::cout << "pi = " << test.probabilityOfSuperiority() << std::endl;
	return 0;
}
