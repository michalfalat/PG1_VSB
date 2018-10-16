#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "vertex.h"

class Surface; // dopøedná deklarace tøídy

/*! \class Triangle
\brief A class representing single triangle in 3D.

\author Tomáš Fabián
\version 1.1
\date 2013-2018
*/
class Triangle
{
public:	
	//! Vıchozí konstruktor.
	/*!
	Inicializuje všechny sloky trojúhelníku na hodnotu nula.
	*/
	Triangle() { }

	//! Obecnı konstruktor.
	/*!
	Inicializuje trojúhelník podle zadanıch hodnot parametrù.

	\param v0 první vrchol trojúhelníka.
	\param v1 druhı vrchol trojúhelníka.
	\param v2 tøetí vrchol trojúhelníka.
	\param surface ukazatel na plochu, jí je trojúhelník èlenem.
	*/
	Triangle( const Vertex & v0, const Vertex & v1, const Vertex & v2, Surface * surface = NULL );

	//! I-tı vrchol trojúhelníka.
	/*!
	\param i index vrcholu trojúhelníka.

	\return I-tı vrchol trojúhelníka.
	*/
	Vertex vertex( const int i );
	
	//! Ukazatel na sí, jí je trojúhelník èlenem.
	/*!
	\return Ukazatel na sí.
	*/	
	Surface * surface();

	//! Vypoète normálu trojúhelníka v zadaném bodì.
	/*!
	 * \param u baricentrická souøadnice bodu, pro kterı má bıt normála interpolována.
	 * \param v baricentrická souøadnice bodu, pro kterı má bıt normála interpolována.
	 *
	 * \return Normála trojúhelníka v zadaném bodì.
	*/
	Vector3 normal(const float u, const float v);

protected:

private:
	Vertex vertices_[3]; /*!< Vrcholy trojúhelníka. Nic jiného tu nesmí bıt, jinak padne VBO v OpenGL! */	
};

#endif
