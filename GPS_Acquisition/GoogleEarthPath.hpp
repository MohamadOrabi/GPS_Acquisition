#ifndef GOOGLEEARTHPATH_HPP_INCLUDED
#define GOOGLEEARTHPATH_HPP_INCLUDED

#include <iostream>	
#include <fstream>	
#include <iomanip>  
#include <unordered_map>

using namespace std;

class GoogleEarthPath
{
public:
	inline			GoogleEarthPath(string, string);	
	inline			~GoogleEarthPath();

	inline void		addPoint(double, double);
private:
	fstream			fileDescriptor;					// File descriptor
};



inline GoogleEarthPath::GoogleEarthPath(string file, string pathName)
{
	fileDescriptor.open(file, ios::out | ios::trunc);	// Delete previous file if it exists	
	if(!fileDescriptor)
		std::cerr << "GoogleEarthPath::\tCould not open file file! " << file << std::endl; 

	fileDescriptor 	<< "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" 	<< "\n" 
			<< "<kml xmlns=\"http://www.opengis.net/kml/2.2\">"	<< "\n" 
			<< "<Document>"						<< "\n" 

			<< "<Style id=\"PathStyle\">" << "\n"
			<< "<LineStyle>" << "\n"
			<< "<color>ff00ff00</color>"
			<< "<width>5</width>" << "\n"
			<< "</LineStyle>" << "\n"
			<< "</Style>" << "\n"

			<< "<Placemark>"					<< "\n" 
			<< "<name>"<< pathName <<"</name>"			<< "\n" 
			<< "<styleUrl>#PathStyle</styleUrl>"	<< "\n"
			<< "<LineString>"					<< "\n" 
			<< "<tessellate>0</tessellate>"				<< "\n" 
			<< "<coordinates>"					<< "\n" ;

}

inline GoogleEarthPath::~GoogleEarthPath()
{
	if(fileDescriptor)
	{
		fileDescriptor	<< "</coordinates>"	<< "\n" 
				<< "</LineString>"	<< "\n" 
				<< "</Placemark>"	<< "\n" 
				<< "</Document>"	<< "\n" 
				<< "</kml>"		<< "\n";
	}
}

inline void GoogleEarthPath::addPoint(double longitude, double latitude)
{
	if(fileDescriptor)
	{
		fileDescriptor 	<< setprecision(13) 
				<< longitude <<"," 
				<< setprecision(13) 
				<< latitude << ",0\n" 
				<< flush;	
	}
}


#endif 


