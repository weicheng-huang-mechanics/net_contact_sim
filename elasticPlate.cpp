#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_radius, 
		double m_Possion, double m_dt, double m_deltaLength)
{
	YoungM = m_YoungM;
	density = m_density;
	radius = m_radius;
	Possion = m_Possion;
	dt = m_dt;

	deltaEdge = m_deltaLength;

	EA = YoungM * M_PI * radius * radius;
	EI = YoungM * M_PI * radius * radius * radius * radius / 4;

	setupGeometry();

	ndof = 3 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(3 * i + 0) = v_nodes[i](0);
		x(3 * i + 1) = v_nodes[i](1);
		x(3 * i + 2) = v_nodes[i](2) + 4.0;
	}
	x0 = x;

	computeEdge();
	computeBending();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }

}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	
	double deltaMass;

	int index1;
	int index2;

	for (int i = 0; i < edgeNum; i++)
	{
		deltaMass = M_PI * radius * radius * density * v_edgeElement[i].refLength / 2;

		index1 = v_edgeElement[i].nv_1;
		index2 = v_edgeElement[i].nv_2;

		massArray(3 * index1 + 0) = massArray(3 * index1 + 0) + deltaMass;
		massArray(3 * index1 + 1) = massArray(3 * index1 + 1) + deltaMass;
		massArray(3 * index1 + 2) = massArray(3 * index1 + 2) + deltaMass;

		massArray(3 * index2 + 0) = massArray(3 * index2 + 0) + deltaMass;
		massArray(3 * index2 + 1) = massArray(3 * index2 + 1) + deltaMass;
		massArray(3 * index2 + 2) = massArray(3 * index2 + 2) + deltaMass;
	}
}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector3d position, int k)
{
	isConstrained[3 * k + 0] = 1;
	isConstrained[3 * k + 1] = 1;
	isConstrained[3 * k + 2] = 1;

	// Store in the constrained dof vector
	x(3 * k + 0) = position(0);
	x(3 * k + 1) = position(1);
	x(3 * k + 2) = position(2);
}

Vector3d elasticPlate::getVertex(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x(3 * i + 0);
	xCurrent(1) = x(3 * i + 1);
	xCurrent(2) = x(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVertexOld(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x0(3 * i + 0);
	xCurrent(1) = x0(3 * i + 1);
	xCurrent(2) = x0(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVelocity(int i)
{
	Vector3d uCurrent;

	uCurrent(0) = ( x(3 * i + 0) - x0(3 * i + 0) ) / dt;
	uCurrent(1) = ( x(3 * i + 1) - x0(3 * i + 1) ) / dt;
	uCurrent(2) = ( x(3 * i + 2) - x0(3 * i + 2) ) / dt;

	return uCurrent;
}

void elasticPlate::updateEdgePair()
{
	for (int i = 0; i < edgeNum; i++)
	{
		v_edgeElement[i].x_1 = getVertex(v_edgeElement[i].nv_1);
		v_edgeElement[i].x_2 = getVertex(v_edgeElement[i].nv_2);
		v_edgeElement[i].edgeLength = (v_edgeElement[i].x_1 - v_edgeElement[i].x_2).norm();
	}
}

void elasticPlate::updateBendingPair()
{
	for (int i = 0; i < bendingNum; i++)
	{
		v_bendingElement[i].x_1 = getVertex(v_bendingElement[i].nv_1);
		v_bendingElement[i].x_2 = getVertex(v_bendingElement[i].nv_2);
		v_bendingElement[i].x_3 = getVertex(v_bendingElement[i].nv_3);

		v_bendingElement[i].e_1 = v_bendingElement[i].x_2 - v_bendingElement[i].x_1;
		v_bendingElement[i].e_2 = v_bendingElement[i].x_3 - v_bendingElement[i].x_2;

		v_bendingElement[i].norm_1 =  v_bendingElement[i].e_1.norm();
		v_bendingElement[i].norm_2 =  v_bendingElement[i].e_2.norm();

		v_bendingElement[i].t_1 = v_bendingElement[i].e_1 / v_bendingElement[i].norm_1;
		v_bendingElement[i].t_2 = v_bendingElement[i].e_2 / v_bendingElement[i].norm_2;
	}
}

void elasticPlate::prepareForIteration()
{
	updateEdgePair();
	updateBendingPair();
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	v_nodes.clear();
    edge.clear();
    bending.clear();
    constraint.clear();

    ifstream inFile1;
	inFile1.open("inputdata/nodesInput.txt");
	double a, b, c;
	nv = 0;
	while(inFile1 >> a >> b >> c)
	{
		Vector3d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;
		xCurrent(2) = c;

		nv = nv + 1;

		v_nodes.push_back(xCurrent);
	}
	inFile1.close();

	ifstream inFile2;
	inFile2.open("inputdata/edgeInput.txt");
	int d, e;
	while(inFile2 >> d >> e)
	{
		Vector2i edgeCurrent;

    	edgeCurrent(0) = d;
    	edgeCurrent(1) = e;

    	edge.push_back(edgeCurrent);
	}
	inFile2.close();

	ifstream inFile3;
	inFile3.open("inputdata/bendingInput.txt");
	int f, g, h;
	while(inFile3 >> f >> g >> h)
	{
		Vector3i bendingCurrent;

    	bendingCurrent(0) = f;
    	bendingCurrent(1) = g; 
    	bendingCurrent(2) = h; 

    	bending.push_back(bendingCurrent);
	}
	inFile3.close();

	ifstream inFile4;
	inFile4.open("inputdata/connerIndex.txt");
	int a1;
	while(inFile4 >> a1)
	{
		int constrainedBC;

		constrainedBC = a1; 

    	constraint.push_back(constrainedBC);
	}
	inFile4.close();



    v_satelliteNodes.clear();
	ifstream inFile7;
	inFile7.open("inputdata/plateNodes.txt");
	double plateNode1, plateNode2, plateNode3;
	while(inFile7 >> plateNode1 >> plateNode2 >> plateNode3)
	{
		Vector3d objectNodes;

		objectNodes(0) = plateNode1; 
		objectNodes(1) = plateNode2; 
		objectNodes(2) = plateNode3; 

    	v_satelliteNodes.push_back(objectNodes);
	}
	inFile7.close();

	v_satelliteMesh.clear();
	ifstream inFile8;
	inFile8.open("inputdata/plateMesh.txt");
	int plateMesh1, plateMesh2, plateMesh3;
	while(inFile8 >> plateMesh1 >> plateMesh2 >> plateMesh3)
	{
		Vector3i objectMesh;

		objectMesh(0) = plateMesh1; 
		objectMesh(1) = plateMesh2;
		objectMesh(2) = plateMesh3; 

    	v_satelliteMesh.push_back(objectMesh);
	}
	inFile8.close();

	v_satelliteLine.clear();
	ifstream inFile9;
	inFile9.open("inputdata/plateLine.txt");
	int plateLine1, plateLine2;
	while(inFile9 >> plateLine1 >> plateLine2)
	{
		Vector2i objectLine;

		objectLine(0) = plateLine1; 
		objectLine(1) = plateLine2;

    	v_satelliteLine.push_back(objectLine);
	}
	inFile9.close();
}

void elasticPlate::computeEdge()
{
	edgeNum = 0;
	v_edgeElement.clear();

	for (int i = 0; i < edge.size(); i++)
	{
		Vector2i edgeCurrent = edge[i];

		edgeElement m_edgeElement;

		m_edgeElement.nv_1 = edgeCurrent(0);
		m_edgeElement.nv_2 = edgeCurrent(1);

		m_edgeElement.x_1 = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2 = getVertex(m_edgeElement.nv_2);

		m_edgeElement.refLength = (m_edgeElement.x_2- m_edgeElement.x_1).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

		m_edgeElement.arrayNum = VectorXi::Zero(6);

		m_edgeElement.arrayNum(0) = 3 * m_edgeElement.nv_1 + 0;
		m_edgeElement.arrayNum(1) = 3 * m_edgeElement.nv_1 + 1;
		m_edgeElement.arrayNum(2) = 3 * m_edgeElement.nv_1 + 2;

		m_edgeElement.arrayNum(3) = 3 * m_edgeElement.nv_2 + 0;
		m_edgeElement.arrayNum(4) = 3 * m_edgeElement.nv_2 + 1;
		m_edgeElement.arrayNum(5) = 3 * m_edgeElement.nv_2 + 2;

		v_edgeElement.push_back(m_edgeElement);

		edgeNum = edgeNum + 1;
	}
}

void elasticPlate::computeBending()
{
	bendingNum = 0;
	v_bendingElement.clear();

	for (int i = 0; i < bending.size(); i++)
	{
		Vector3i bendingCurrent = bending[i];

		bendingElement m_bendingElement;

		m_bendingElement.nv_1 = bendingCurrent(0);
		m_bendingElement.nv_2 = bendingCurrent(1);
		m_bendingElement.nv_3 = bendingCurrent(2);

		m_bendingElement.x_1 = getVertex(m_bendingElement.nv_1);
		m_bendingElement.x_2 = getVertex(m_bendingElement.nv_2);
		m_bendingElement.x_3 = getVertex(m_bendingElement.nv_3);

		m_bendingElement.e_1 = m_bendingElement.x_2 - m_bendingElement.x_1;
		m_bendingElement.e_2 = m_bendingElement.x_3 - m_bendingElement.x_2;

		m_bendingElement.norm_1 =  m_bendingElement.e_1.norm();
		m_bendingElement.norm_2 =  m_bendingElement.e_2.norm();

		m_bendingElement.voroniLength = (m_bendingElement.norm_1 + m_bendingElement.norm_2) / 2;

		m_bendingElement.t_1 = m_bendingElement.e_1 / m_bendingElement.norm_1;
		m_bendingElement.t_2 = m_bendingElement.e_2 / m_bendingElement.norm_2;

		m_bendingElement.nBar = m_bendingElement.t_2 - m_bendingElement.t_1;

		m_bendingElement.arrayNum = VectorXi::Zero(9);

		m_bendingElement.arrayNum(0) = 3 * m_bendingElement.nv_1 + 0;
		m_bendingElement.arrayNum(1) = 3 * m_bendingElement.nv_1 + 1;
		m_bendingElement.arrayNum(2) = 3 * m_bendingElement.nv_1 + 2;

		m_bendingElement.arrayNum(3) = 3 * m_bendingElement.nv_2 + 0;
		m_bendingElement.arrayNum(4) = 3 * m_bendingElement.nv_2 + 1;
		m_bendingElement.arrayNum(5) = 3 * m_bendingElement.nv_2 + 2;

		m_bendingElement.arrayNum(6) = 3 * m_bendingElement.nv_3 + 0;
		m_bendingElement.arrayNum(7) = 3 * m_bendingElement.nv_3 + 1;
		m_bendingElement.arrayNum(8) = 3 * m_bendingElement.nv_3 + 2;

		v_bendingElement.push_back(m_bendingElement);

		bendingNum = bendingNum + 1;
	}
}