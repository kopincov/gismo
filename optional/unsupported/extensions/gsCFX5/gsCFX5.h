/** @file gsCFX5.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moeller
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include <cfxImport.h>
#include <cfxExport.h>

#include <gsCore/gsDebug.h>
#include <gsUtils/gsUtils.h>

namespace gismo {

  // Forward declarations
  class gsCFX5;
  class gsCFX5Boundary;
  class gsCFX5Domain;
  class gsCFX5Element;
  class gsCFX5ElementInVolumeConstIterator;
  class gsCFX5ElementInVolumeConstReverseIterator;
  class gsCFX5ElementInVolumeIterator;
  class gsCFX5ElementInVolumeReverseIterator;
  class gsCFX5Face;
  class gsCFX5FaceInRegionConstIterator;
  class gsCFX5FaceInRegionConstReverseIterator;
  class gsCFX5Node;
  class gsCFX5NodeInRegionConstIterator;
  class gsCFX5NodeInRegionConstReverseIterator;
  class gsCFX5NodeInRegionIterator;
  class gsCFX5NodeInRegionReverseIterator;
  class gsCFX5NodeInVolumeConstIterator;
  class gsCFX5NodeInVolumeConstReverseIterator;
  class gsCFX5NodeInVolumeIterator;
  class gsCFX5NodeInVolumeReverseIterator;
  class gsCFX5Region;
  class gsCFX5Variable;
  class gsCFX5Volume;
  template<typename T> class gsCFX5GenericConstIterator;
  template<typename T> class gsCFX5GenericConstReverseIterator;
  template<typename T> class gsCFX5GenericIterator;
  template<typename T> class gsCFX5GenericReverseIterator;

  /**
     \brief Deep-copy of CFX5 pointer

     The CFX5 mesh import/export API only returns pointers to dynamic
     data structures which are stored internally and are deleted once
     an import/export file is closed. As a remedy, a deep-copy of the
     dynamic data structures must be performed in order to be able to
     import/export several files at a time.
  */
  template<typename S, typename T=S>
  T* cfxCopy(const S* data, size_t size)
  {
#if !defined(NDEBUG)
    gsDebug << "cfxCopy type " << typeid(S).name()
            <<" of size " << size
            << " to type " << typeid(T).name() << std::endl;
#endif
    T* t = new T[size];
    for (size_t i=0; i<size; i++)
      t[i] = (T)data[i];

    return t;
  }

  /**
     \brief A CFX5 generic iterator.
  */
  template<typename T>
  class gsCFX5GenericIterator
  {
  public:
    // Iterator traits, previously from std::iterator
    typedef gsCFX5GenericIterator     self_type;
    typedef T                         value_type;
    typedef T&                        reference;
    typedef T*                        pointer;
    typedef std::forward_iterator_tag iterator_category;
    typedef ptrdiff_t                 difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5GenericIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5GenericIterator(pointer ptr)
      :
      ptr(ptr)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return *ptr; }

    /// \brief Pointer operator
    value_type* operator->() { return ptr; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) { return ptr == rhs.ptr; }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) { return ptr != rhs.ptr; }

  private:
    pointer ptr = nullptr;
  };

  /**
     \brief A CFX5 generic const-iterator.
  */
  template<typename T>
  class gsCFX5GenericConstIterator
  {
  public:
    // Iterator traits, previously from std::const_iterator
    typedef gsCFX5GenericConstIterator self_type;
    typedef T                          value_type;
    typedef T&                         reference;
    typedef T*                         pointer;
    typedef std::forward_iterator_tag  iterator_category;
    typedef ptrdiff_t                  difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5GenericConstIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5GenericConstIterator(pointer ptr)
      :
      ptr(ptr)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return *ptr; }

    /// \brief Pointer operator
    const value_type* operator->() const { return ptr; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const { return ptr == rhs.ptr; }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const { return ptr != rhs.ptr; }

  private:
    pointer ptr = nullptr;
  };

  /**
     \brief A CFX5 generic reverse-iterator.
  */
  template<typename T>
  class gsCFX5GenericReverseIterator
  {
  public:
    // Iterator traits, previously from std::reverse_iterator
    typedef gsCFX5GenericReverseIterator self_type;
    typedef T                            value_type;
    typedef T&                           reference;
    typedef T*                           pointer;
    typedef std::forward_iterator_tag    iterator_category;
    typedef ptrdiff_t                    difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5GenericReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5GenericReverseIterator(pointer ptr)
      :
      ptr(ptr)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return *ptr; }

    /// \brief Pointer operator
    value_type* operator->() { return ptr; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) { return ptr == rhs.ptr; }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) { return ptr != rhs.ptr; }

  private:
    pointer ptr = nullptr;
  };

  /**
     \brief A CFX5 generic const-reverse-iterator.
  */
  template<typename T>
  class gsCFX5GenericConstReverseIterator
  {
  public:
    // Iterator traits, previously from std::const_reverse_iterator
    typedef gsCFX5GenericConstReverseIterator self_type;
    typedef T                                 value_type;
    typedef T&                                reference;
    typedef T*                                pointer;
    typedef std::forward_iterator_tag         iterator_category;
    typedef ptrdiff_t                         difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5GenericConstReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5GenericConstReverseIterator(pointer ptr)
      :
      ptr(ptr)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return *ptr; }

    /// \brief Pointer operator
    const value_type* operator->() const { return ptr; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const { return ptr == rhs.ptr; }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const { return ptr != rhs.ptr; }

  private:
    pointer ptr = nullptr;
  };

  /**
     \brief CFX5 element types
  */
  enum class gsCFX5ElementType : int
  {
    /* Documentation from cfxExport.h

       cfxELEM_TET     4   tet element     (4 nodes )
       cfxELEM_PYR     5   pyramid element (5 nodes )
       cfxELEM_WDG     6   wedge element   (6 nodes)
       cfxELEM_HEX     8   hex element     (8 nodes)
    */

      Tetrahedral = cfxELEM_TET,
      Pyramid     = cfxELEM_PYR,
      Prism       = cfxELEM_WDG,
      Hexahedral  = cfxELEM_HEX
      };

  /**
     \brief A CFX5 element (wrapper for cfxElement).
  */
  class gsCFX5Element : public cfxElement
  {
    /* Documentation from cfxExport.h

       typedef struct cfxElement {
       int type;
       int *nodeid;
       } cfxElement;
    */

  public:
    /// \brief Constructor (default)
    gsCFX5Element()
    {
      type   = 0;
      nodeid = nullptr;
    }

    /// \brief Constructor (copy)
    gsCFX5Element(const gsCFX5Element& other)
    {
#if !defined(NDEBUG)
      gsDebug << "Copying gsCFX5Element\n";
#endif
      type = (int)other.getElementType();
      nodeid = new int[type];
      for (int i=0; i<type; i++)
        nodeid[i] = other.getNodeID(i);
    }

    /// \brief Constructor (move)
    gsCFX5Element(gsCFX5Element&& other)
    {
#if !defined(NDEBUG)
      gsDebug << "Moving gsCFX5Element\n";
#endif
      type = (int)other.getElementType();
      nodeid = other.getNodeIDs();;
      other.nodeid = nullptr;
    }

    /// \brief Destructor
    ~gsCFX5Element()
    {
      type = 0;
      delete[] nodeid;
    }

    /// \brief Copy assignment operator
    gsCFX5Element& operator=(const gsCFX5Element& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Copy gsCFX5Element\n";
#endif
          type = (int)other.getElementType();
          delete[] nodeid;
          nodeid = new int[type];
          for (int i=0; i<type; i++)
            nodeid[i] = other.getNodeID(i);
        }
      return *this;
    }

    /// \brief Move assignment operator
    gsCFX5Element& operator=(gsCFX5Element&& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Move assigning gsCFX5Element\n";
#endif
          type = (int)other.getElementType();
          delete[] nodeid;
          nodeid = other.getNodeIDs();
          other.nodeid = nullptr;
        }
      return *this;
    }

    /// \brief Returns swapped element
    gsCFX5Element& swap(gsCFX5Element& other)
    {
      std::swap(type,   other.type);
      std::swap(nodeid, other.nodeid);
      return *this;
    }

  public:
    /// \brief Returns the type of the element
    const gsCFX5ElementType getElementType() const
    { return (gsCFX5ElementType)type; }

    /// \brief Sets the type of the element
    void setElementType(gsCFX5ElementType type)
    {
      this->type = (int)type;
      delete[] nodeid;
      nodeid = new int[this->type];
    }

    /// \brief Returns constant pointer to the node IDs of the element
    const int* getNodeIDs() const
    { return nodeid; }

    /// \brief Returns pointer to the node IDs of the element
    int* getNodeIDs()
    { return nodeid; }

    /// \brief Returns constant reference to the idx-th node ID of the element
    const int& getNodeID(const int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < (int)getElementType()), "Invalid node index.");
      return nodeid[idx];
    }

    /// \brief Returns reference to the idx-th node ID of the element
    int& getNodeID(const int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < (int)getElementType()), "Invalid node index.");
      return nodeid[idx];
    }

    /// \brief Compares two elements for equality
    bool operator==(gsCFX5Element& other)
    {
      if ((gsCFX5ElementType)type != other.getElementType())
        return false;

      for (int idx=0; idx<(int)type; idx++)
        if (nodeid[idx] != other.getNodeID(idx))
          return false;

      return true;
    }

    /// \brief Compares two elements for inequality
    bool operator!=(gsCFX5Element& other)
    {
      if ((gsCFX5ElementType)type != other.getElementType())
        return true;

      for (int idx=0; idx<(int)type; idx++)
        if (nodeid[idx] != other.getNodeID(idx))
          return true;

      return false;
    }

    /**
       \brief Returns elements with node ordering swapped

       The CFX5 mesh import/export API uses different formats for
       storing the nodes in hexahedral elements. The internal storage
       of nodes in the gsCFX5Element class adopts the cfxExport
       format. This method converts between the two formats. By
       calling it twice, the original format is recovered.
    */
    gsCFX5Element swapExportImport() const
    {
      gsCFX5Element element(*this);

      // Swap nodes
      if ((gsCFX5ElementType)type == gsCFX5ElementType::Hexahedral)
        {
          element.getNodeID(2) = getNodeID(3);
          element.getNodeID(3) = getNodeID(2);

          element.getNodeID(6) = getNodeID(7);
          element.getNodeID(7) = getNodeID(6);
        }

      return element;
    }

    /// \brief Returns hash of the element
    size_t getHash() const
    {
      switch (getElementType())
        {
        case gsCFX5ElementType::Tetrahedral :
          return std::hash<std::string>{}(util::to_string(getNodeID(0))+
                                          util::to_string(getNodeID(1))+
                                          util::to_string(getNodeID(2))+
                                          util::to_string(getNodeID(3)));

        case gsCFX5ElementType::Pyramid :
          return std::hash<std::string>{}(util::to_string(getNodeID(0))+
                                          util::to_string(getNodeID(1))+
                                          util::to_string(getNodeID(2))+
                                          util::to_string(getNodeID(3))+
                                          util::to_string(getNodeID(4)));

        case gsCFX5ElementType::Prism :
          return std::hash<std::string>{}(util::to_string(getNodeID(0))+
                                          util::to_string(getNodeID(1))+
                                          util::to_string(getNodeID(2))+
                                          util::to_string(getNodeID(3))+
                                          util::to_string(getNodeID(4))+
                                          util::to_string(getNodeID(5)));

        case gsCFX5ElementType::Hexahedral :
          return std::hash<std::string>{}(util::to_string(getNodeID(0))+
                                          util::to_string(getNodeID(1))+
                                          util::to_string(getNodeID(2))+
                                          util::to_string(getNodeID(3))+
                                          util::to_string(getNodeID(4))+
                                          util::to_string(getNodeID(5))+
                                          util::to_string(getNodeID(6))+
                                          util::to_string(getNodeID(7)));

        default:
          GISMO_ERROR("Invalid element type!");
        }
    }

    /// \brief Prints the element
    std::ostream& print(std::ostream& os) const
    {
      switch (getElementType())
        {
        case gsCFX5ElementType::Tetrahedral :
          os << "TET ("
             << getNodeID(0) << ","
             << getNodeID(1) << ","
             << getNodeID(2) << ","
             << getNodeID(3) << ")";
          break;

        case gsCFX5ElementType::Pyramid :
          os << "PYR ("
             << getNodeID(0) << ","
             << getNodeID(1) << ","
             << getNodeID(2) << ","
             << getNodeID(3) << ","
             << getNodeID(4) << ")";
          break;

        case gsCFX5ElementType::Prism :
          os << "PRI ("
             << getNodeID(0) << ","
             << getNodeID(1) << ","
             << getNodeID(2) << ","
             << getNodeID(3) << ","
             << getNodeID(4) << ","
             << getNodeID(5) << ")";
          break;

        case gsCFX5ElementType::Hexahedral :
          os << "HEX ("
             << getNodeID(0) << ","
             << getNodeID(1) << ","
             << getNodeID(2) << ","
             << getNodeID(3) << ","
             << getNodeID(4) << ","
             << getNodeID(5) << ","
             << getNodeID(6) << ","
             << getNodeID(7) << ")";
          break;

        default:
          GISMO_ERROR("Invalid element type!");
        }
      return os;
    }

    /// \brief Prints the element as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);
      return os.str();
    }

    /// \brief Prints the difference between the element and the other
    /// element as a string
    std::string diff(const gsCFX5Element& other) const
    {
      std::ostringstream os, str, other_str;
      print(str);
      other.print(other_str);

      if (str.str() != other_str.str())
        os << "< " << str.str() << "\n---\n> " << other_str.str() << "\n";

      return os.str();
    }
  };

  /// \brief Output CFX5 element object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Element& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 node (wrapper for cfxNode).
  */
  class gsCFX5Node : public cfxNode
  {
    /* Documentation from cfxExport.h

       typedef struct cfxNode {
       double x, y, z;
       } cfxNode;
    */

  public:
    /// \brief Returns constant reference to X-coordinate of the node
    const double& getX() const
    { return x; }

    /// \brief Returns constant reference to Y-coordinate of the node
    const double& getY() const
    { return y; }

    /// \brief Returns constant reference to Z-coordinate of the node
    const double& getZ() const
    { return z; }

    /// \brief Returns reference to X-coordinate of the node
    double& getX()
    { return x; }

    /// \brief Returns reference to Y-coordinate of the node
    double& getY()
    { return y; }

    /// \brief Returns reference to Z-coordinate of the node
    double& getZ()
    { return z; }

    /// \brief Compares two nodes for equality
    bool operator==(const gsCFX5Node& other)
    {
      return ((x==other.getX()) && (y==other.getY()) && (z==other.getZ()));
    }

    /// \brief Compares two nodes for inequality
    bool operator!=(const gsCFX5Node& other)
    {
      return ((x!=other.getX()) || (y!=other.getY()) || (z!=other.getZ()));
    }

    /// \brief Returns hash of the node
    size_t getHash() const
    {
      return std::hash<std::string>{}(util::to_string(x)+
                                      util::to_string(y)+
                                      util::to_string(z));
    }

    /// \brief Prints the node
    std::ostream& print(std::ostream& os) const
    {
      os << "NODE ("
         << getX() << ","
         << getY() << ","
         << getZ() << ")";
      return os;
    }

    /// \brief Prints the node as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);
      return os.str();
    }

    /// \brief Prints the difference between the node and the other
    /// node as a string
    std::string diff(const gsCFX5Node& other) const
    {
      std::ostringstream os, str, other_str;
      print(str);
      other.print(other_str);

      if (str.str() != other_str.str())
        os << "< " << str.str() << "\n---\n> " << other_str.str() << "\n";

      return os.str();
    }
  };

  /// \brief Output CFX5 node object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Node& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 node-in-region iterator.
  */
  class gsCFX5NodeInRegionIterator
  {
  public:
    // Iterator traits, previously from std::iterator
    typedef gsCFX5NodeInRegionIterator self_type;
    typedef gsCFX5Node                 value_type;
    typedef gsCFX5Node&                reference;
    typedef gsCFX5Node*                pointer;
    typedef std::forward_iterator_tag  iterator_category;
    typedef ptrdiff_t                  difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInRegionIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInRegionIterator(int* ptr, gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    value_type* operator->() { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs)
    { return (ptr == rhs.ptr) && (ptr_node == rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs)
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-region const-iterator.
  */
  class gsCFX5NodeInRegionConstIterator
  {
  public:
    // Iterator traits, previously from std::const_iterator
    typedef gsCFX5NodeInRegionConstIterator self_type;
    typedef gsCFX5Node                      value_type;
    typedef gsCFX5Node&                     reference;
    typedef gsCFX5Node*                     pointer;
    typedef std::forward_iterator_tag       iterator_category;
    typedef ptrdiff_t                       difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInRegionConstIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInRegionConstIterator(int* ptr, const gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    const value_type* operator->() const { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_node == rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    const gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-region reverse-iterator.
  */
  class gsCFX5NodeInRegionReverseIterator
  {
  public:
    // Iterator traits, previously from std::reverse_iterator
    typedef gsCFX5NodeInRegionReverseIterator self_type;
    typedef gsCFX5Node                        value_type;
    typedef gsCFX5Node&                       reference;
    typedef gsCFX5Node*                       pointer;
    typedef std::forward_iterator_tag         iterator_category;
    typedef ptrdiff_t                         difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInRegionReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInRegionReverseIterator(int* ptr, gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    value_type* operator->() { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-region const-reverse-iterator.
  */
  class gsCFX5NodeInRegionConstReverseIterator
  {
  public:
    // Iterator traits, previously from std::const_reverse_iterator
    typedef gsCFX5NodeInRegionConstReverseIterator self_type;
    typedef gsCFX5Node                             value_type;
    typedef gsCFX5Node&                            reference;
    typedef gsCFX5Node*                            pointer;
    typedef std::forward_iterator_tag              iterator_category;
    typedef ptrdiff_t                              difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInRegionConstReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInRegionConstReverseIterator(int* ptr, const gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    const value_type* operator->() const { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    const gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-volume iterator.
  */
  class gsCFX5NodeInVolumeIterator
  {
  public:
    // Iterator traits, previously from std::iterator
    typedef gsCFX5NodeInVolumeIterator self_type;
    typedef gsCFX5Node                 value_type;
    typedef gsCFX5Node&                reference;
    typedef gsCFX5Node*                pointer;
    typedef std::forward_iterator_tag  iterator_category;
    typedef ptrdiff_t                  difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInVolumeIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInVolumeIterator(int* ptr, gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    value_type* operator->() { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs)
    { return (ptr == rhs.ptr) && (ptr_node == rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs)
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-volume const-iterator.
  */
  class gsCFX5NodeInVolumeConstIterator
  {
  public:
    // Iterator traits, previously from std::const_iterator
    typedef gsCFX5NodeInVolumeConstIterator self_type;
    typedef gsCFX5Node                      value_type;
    typedef gsCFX5Node&                     reference;
    typedef gsCFX5Node*                     pointer;
    typedef std::forward_iterator_tag       iterator_category;
    typedef ptrdiff_t                       difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInVolumeConstIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInVolumeConstIterator(int* ptr, const gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    const value_type* operator->() const { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_node == rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    const gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-volume reverse-iterator.
  */
  class gsCFX5NodeInVolumeReverseIterator
  {
  public:
    // Iterator traits, previously from std::reverse_iterator
    typedef gsCFX5NodeInVolumeReverseIterator self_type;
    typedef gsCFX5Node                        value_type;
    typedef gsCFX5Node&                       reference;
    typedef gsCFX5Node*                       pointer;
    typedef std::forward_iterator_tag         iterator_category;
    typedef ptrdiff_t                         difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInVolumeReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInVolumeReverseIterator(int* ptr, gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    value_type* operator->() { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 node-in-volume const-reverse-iterator.
  */
  class gsCFX5NodeInVolumeConstReverseIterator
  {
  public:
    // Iterator traits, previously from std::const_reverse_iterator
    typedef gsCFX5NodeInVolumeConstReverseIterator self_type;
    typedef gsCFX5Node                             value_type;
    typedef gsCFX5Node&                            reference;
    typedef gsCFX5Node*                            pointer;
    typedef std::forward_iterator_tag              iterator_category;
    typedef ptrdiff_t                              difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5NodeInVolumeConstReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5NodeInVolumeConstReverseIterator(int* ptr, const gsCFX5Node* ptr_node)
      :
      ptr(ptr),
      ptr_node(ptr_node)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_node = other.ptr_node; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return ptr_node[*ptr-1]; }

    /// \brief Pointer operator
    const value_type* operator->() const { return &ptr_node[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_node != rhs.ptr_node); }

  private:
    int* ptr = nullptr;
    const gsCFX5Node* ptr_node = nullptr;
  };

  /**
     \brief A CFX5 element-in-volume iterator.
  */
  class gsCFX5ElementInVolumeIterator
  {
  public:
    // Iterator traits, previously from std::iterator
    typedef gsCFX5ElementInVolumeIterator self_type;
    typedef gsCFX5Element                 value_type;
    typedef gsCFX5Element&                reference;
    typedef gsCFX5Element*                pointer;
    typedef std::forward_iterator_tag     iterator_category;
    typedef ptrdiff_t                     difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5ElementInVolumeIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5ElementInVolumeIterator(int* ptr, gsCFX5Element* ptr_element)
      :
      ptr(ptr),
      ptr_element(ptr_element)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_element = other.ptr_element; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return ptr_element[*ptr-1]; }

    /// \brief Pointer operator
    value_type* operator->() { return &ptr_element[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs)
    { return (ptr == rhs.ptr) && (ptr_element == rhs.ptr_element); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs)
    { return (ptr != rhs.ptr) || (ptr_element != rhs.ptr_element); }

  private:
    int* ptr = nullptr;
    gsCFX5Element* ptr_element = nullptr;
  };

  /**
     \brief A CFX5 element-in-volume const-iterator.
  */
  class gsCFX5ElementInVolumeConstIterator
  {
  public:
    // Iterator traits, previously from std::const_iterator
    typedef gsCFX5ElementInVolumeConstIterator self_type;
    typedef gsCFX5Element                      value_type;
    typedef gsCFX5Element&                     reference;
    typedef gsCFX5Element*                     pointer;
    typedef std::forward_iterator_tag          iterator_category;
    typedef ptrdiff_t                          difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5ElementInVolumeConstIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5ElementInVolumeConstIterator(int* ptr, const gsCFX5Element* ptr_element)
      :
      ptr(ptr),
      ptr_element(ptr_element)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_element = other.ptr_element; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return ptr_element[*ptr-1]; }

    /// \brief Pointer operator
    const value_type* operator->() const { return &ptr_element[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_element == rhs.ptr_element); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_element != rhs.ptr_element); }

  private:
    int* ptr = nullptr;
    const gsCFX5Element* ptr_element = nullptr;
  };

  /**
     \brief A CFX5 element-in-volume reverse-iterator.
  */
  class gsCFX5ElementInVolumeReverseIterator
  {
  public:
    // Iterator traits, previously from std::reverse_iterator
    typedef gsCFX5ElementInVolumeReverseIterator self_type;
    typedef gsCFX5Element                        value_type;
    typedef gsCFX5Element&                       reference;
    typedef gsCFX5Element*                       pointer;
    typedef std::forward_iterator_tag            iterator_category;
    typedef ptrdiff_t                            difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5ElementInVolumeReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5ElementInVolumeReverseIterator(int* ptr, gsCFX5Element* ptr_element)
      :
      ptr(ptr),
      ptr_element(ptr_element)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_element = other.ptr_element; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    value_type& operator*() { return ptr_element[*ptr-1]; }

    /// \brief Pointer operator
    value_type* operator->() { return &ptr_element[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_element == rhs.ptr_element); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_element != rhs.ptr_element); }

  private:
    int* ptr = nullptr;
    gsCFX5Element* ptr_element = nullptr;
  };

  /**
     \brief A CFX5 element-in-volume const-reverse-iterator.
  */
  class gsCFX5ElementInVolumeConstReverseIterator
  {
  public:
    // Iterator traits, previously from std::const_reverse_iterator
    typedef gsCFX5ElementInVolumeConstReverseIterator self_type;
    typedef gsCFX5Element                             value_type;
    typedef gsCFX5Element&                            reference;
    typedef gsCFX5Element*                            pointer;
    typedef std::forward_iterator_tag                 iterator_category;
    typedef ptrdiff_t                                 difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5ElementInVolumeConstReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5ElementInVolumeConstReverseIterator(int* ptr, const gsCFX5Element* ptr_element)
      :
      ptr(ptr),
      ptr_element(ptr_element)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_element = other.ptr_element; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    const value_type& operator*() const { return ptr_element[*ptr-1]; }

    /// \brief Pointer operator
    const value_type* operator->() const { return &ptr_element[*ptr-1]; }

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_element == rhs.ptr_element); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_element != rhs.ptr_element); }

  private:
    int* ptr = nullptr;
    const gsCFX5Element* ptr_element = nullptr;
  };

  /**
     \brief CFX5 face types
  */
  enum class gsCFX5FaceType : int
  {
    Triangle      = 3,
      Quadrilateral = 4
      };

  /**
     \brief A CFX5 face.
  */
  class gsCFX5Face
  {
  public:
    /// \brief Constructor (default)
    gsCFX5Face()
      :
      type{0},
      nodeid{nullptr}
    {
    }

    /// \brief Constructor
    gsCFX5Face(int type, const int* nodeid)
      :
      type{type},
      nodeid{new int[type]}
    {
      for (int i=0; i<type; i++)
        this->nodeid[i] = nodeid[i];
    }

    /// \brief Constructor
    gsCFX5Face(gsCFX5FaceType type, const int* nodeid)
      : gsCFX5Face((int)type, nodeid)
    {
    }

    /// \brief Constructor
    gsCFX5Face(int node0, int node1, int node2)
      :
      type{(int)gsCFX5FaceType::Triangle},
      nodeid{new int[3]}
    {
      nodeid[0] = node0;
      nodeid[1] = node1;
      nodeid[2] = node2;
    }

    /// \brief Constructor
    gsCFX5Face(int node0, int node1, int node2, int node3)
      :
      type{(int)gsCFX5FaceType::Quadrilateral},
      nodeid{new int[4]}
    {
      nodeid[0] = node0;
      nodeid[1] = node1;
      nodeid[2] = node2;
      nodeid[3] = node3;
    }

    /// \brief Constructor (copy)
    gsCFX5Face(const gsCFX5Face& other)
      :
      gsCFX5Face(other.getFaceType(), other.getNodeIDs())
    {
#if !defined(NDEBUG)
      gsDebug << "Copying gsCFX5Face\n";
#endif
    }

    /// \brief Constructor (move)
    gsCFX5Face(gsCFX5Face&& other)
      :
      type{(int)other.getFaceType()}
    {
#if !defined(NDEBUG)
      gsDebug << "Moving gsCFX5Face\n";
#endif
      nodeid = other.getNodeIDs();
      other.nodeid = nullptr;
    }

    /// \brief Destructor
    ~gsCFX5Face()
    {
      type = 0;
      delete[] nodeid;
      nodeid = nullptr;
    }

    /// \brief Copy assignment operator
    gsCFX5Face& operator=(const gsCFX5Face& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Copy assigning gsCFX5Face\n";
#endif
          type = (int)other.getFaceType();
          delete[] nodeid;
          nodeid = new int[type];
          for (int i=0; i<type; i++)
            this->nodeid[i] = other.getNodeID(i);
        }
      return *this;
    }

    /// \brief Move assignment operator
    gsCFX5Face& operator=(gsCFX5Face&& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Move assigning gsCFX5Face\n";
#endif
          type = (int)other.getFaceType();
          delete[] nodeid;
          nodeid = other.getNodeIDs();
          other.nodeid = nullptr;
        }
      return *this;
    }

    /// Returns swapped face
    gsCFX5Face& swap(gsCFX5Face& other)
    {
      std::swap(type,   other.type);
      std::swap(nodeid, other.nodeid);
      return *this;
    }

  public:
    /// \brief Returns the type of the face
    const gsCFX5FaceType getFaceType() const
    { return (gsCFX5FaceType)type; }

    /// \brief Returns constant pointer to the node IDs of the face
    const int* getNodeIDs() const
    { return nodeid; }

    /// \brief Returns pointer to the node IDs of the face
    int* getNodeIDs()
    { return nodeid; }

    /// \brief Returns constant reference to the idx-th node ID of the face
    const int& getNodeID(const int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < (int)getFaceType()), "Invalid node index.");
      return nodeid[idx];
    }

    /// \brief Returns reference to the idx-th node ID of the face
    int& getNodeID(const int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < (int)getFaceType()), "Invalid node index.");
      return nodeid[idx];
    }

    /// \brief Returns hash of the face
    size_t getHash() const
    {
      switch (getFaceType())
        {
        case gsCFX5FaceType::Triangle :
          return std::hash<std::string>{}(util::to_string(getNodeID(0))+
                                          util::to_string(getNodeID(1))+
                                          util::to_string(getNodeID(2)));

        case gsCFX5FaceType::Quadrilateral :
          return std::hash<std::string>{}(util::to_string(getNodeID(0))+
                                          util::to_string(getNodeID(1))+
                                          util::to_string(getNodeID(2))+
                                          util::to_string(getNodeID(3)));

        default:
          GISMO_ERROR("Invalid element type!");
        }
    }

    /// \brief Prints the face
    std::ostream& print(std::ostream& os) const
    {
      switch (getFaceType())
        {
        case gsCFX5FaceType::Triangle :
          os << "TRI ("
             << getNodeID(0) << ","
             << getNodeID(1) << ","
             << getNodeID(2) << ")";
          break;

        case gsCFX5FaceType::Quadrilateral :
          os << "QUA ("
             << getNodeID(0) << ","
             << getNodeID(1) << ","
             << getNodeID(2) << ","
             << getNodeID(3) << ")";
          break;

        default:
          GISMO_ERROR("Invalid face type!");
        }
      return os;
    }

    /// \brief Prints the face as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);
      return os.str();
    }

    /// \brief Prints the difference between the face and the other
    /// face as a string
    std::string diff(const gsCFX5Face& other) const
    {
      std::ostringstream os, str, other_str;
      print(str);
      other.print(other_str);

      if (str.str() != other_str.str())
        os << "< " << str.str() << "\n---\n> " << other_str.str() << "\n";

      return os.str();
    }

  private:
    // Number of nodes in the face
    int type = 0;

    // List of node IDs in the face
    int* nodeid = nullptr;
  };

  /// \brief Output CFX5 face object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Face& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 face-in-region const-iterator.
  */
  class gsCFX5FaceInRegionConstIterator
  {
  public:
    // Iterator traits, previously from std::const_iterator
    typedef gsCFX5FaceInRegionConstIterator self_type;
    typedef gsCFX5Face                      value_type;
    typedef gsCFX5Face&                     reference;
    typedef gsCFX5Face*                     pointer;
    typedef std::forward_iterator_tag       iterator_category;
    typedef ptrdiff_t                       difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5FaceInRegionConstIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5FaceInRegionConstIterator(int ptr, const gsCFX5Region* ptr_region)
      :
      ptr(ptr),
      ptr_region(ptr_region)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_region = other.ptr_region; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr++; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr++; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr--; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr--; return i; }

    /// \brief Dereference operator
    const value_type operator*() const;

    /// \brief Pointer operator
    const value_type operator->() const;

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_region == rhs.ptr_region); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_region != rhs.ptr_region); }

  private:
    int ptr = 0;
    const gsCFX5Region* ptr_region = nullptr;
  };

  /**
     \brief A CFX5 face-in-region const-reverse-iterator.
  */
  class gsCFX5FaceInRegionConstReverseIterator
  {
  public:
    // Iterator traits, previously from std::const_reverse_iterator
    typedef gsCFX5FaceInRegionConstReverseIterator self_type;
    typedef gsCFX5Face                             value_type;
    typedef gsCFX5Face&                            reference;
    typedef gsCFX5Face*                            pointer;
    typedef std::forward_iterator_tag              iterator_category;
    typedef ptrdiff_t                              difference_type;

  public:
    /// \brief Constructor (default)
    gsCFX5FaceInRegionConstReverseIterator() = delete;

    /// \brief Constructor
    explicit gsCFX5FaceInRegionConstReverseIterator(int ptr, const gsCFX5Region* ptr_region)
      :
      ptr(ptr),
      ptr_region(ptr_region)
    {}

    /// \brief Copy assignment operator
    self_type operator=(const self_type& other)
    { ptr = other.ptr; ptr_region = other.ptr_region; return *this; }

    /// \brief Pre-increment operator
    self_type operator++() { ptr--; return *this; }

    /// \brief Post-increment operator
    self_type operator++(int) { self_type i = *this; ptr--; return i; }

    /// \brief Pre-decrement operator
    self_type operator--() { ptr++; return *this; }

    /// \brief Post-decrement operator
    self_type operator--(int) { self_type i = *this; ptr++; return i; }

    /// \brief Dereference operator
    const value_type operator*() const;

    /// \brief Comparison of equality
    bool operator==(const self_type& rhs) const
    { return (ptr == rhs.ptr) && (ptr_region == rhs.ptr_region); }

    /// \brief Comparison of inequality
    bool operator!=(const self_type& rhs) const
    { return (ptr != rhs.ptr) || (ptr_region != rhs.ptr_region); }

  private:
    int ptr = 0;
    const gsCFX5Region* ptr_region = nullptr;
  };

  /**
     \brief A CFX5 region.
  */
  class gsCFX5Region
  {
  public:
    /// \brief Constructor (default)
    gsCFX5Region()
      :
      id{0},
      name{"Default Region"},
      nodes{0},
      faces{0},
      m_node_id{nullptr},
      m_face_id{nullptr},
      m_domain{nullptr}
    {}

    /// \brief Constructor (with pointer to underlying domain)
    gsCFX5Region(gsCFX5Domain* m_domain)
      :
      id{0},
      name{"Default Region"},
      nodes{0},
      faces{0},
      m_node_id{nullptr},
      m_face_id{nullptr},
      m_domain{m_domain}
    {}

    /// \brief Constructor (copy)
    gsCFX5Region(const gsCFX5Region& other)
      :
      id{other.id},
      name{other.name},
      nodes{other.nodes},
      faces{other.faces},
      m_node_id{new int[other.nodes]},
      m_face_id{new int[other.faces]},
      m_domain{other.m_domain}
    {
#if !defined(NDEBUG)
      gsDebug << "Copying gsCFX5Region\n";
#endif
      for (int i=0; i<nodes; i++)
        m_node_id[i] = other.m_node_id[i];

      for (int i=0; i<faces; i++)
        m_face_id[i] = other.m_face_id[i];
    }

    /// \brief Constructor (move)
    gsCFX5Region(gsCFX5Region&& other)
      :
      id{other.id},
      name{other.name},
      nodes{other.nodes},
      faces{other.faces},
      m_node_id{other.m_node_id},
      m_face_id{other.m_face_id},
      m_domain{other.m_domain}
    {
#if !defined(NDEBUG)
      gsDebug << "Moving gsCFX5Region\n";
#endif
      other.id        = 0;
      other.name      = "Default Region";
      other.nodes     = 0;
      other.faces     = 0;
      other.m_node_id = nullptr;
      other.m_face_id = nullptr;
      other.m_domain  = nullptr;
    }

    /// \brief Destructor
    ~gsCFX5Region()
    {
      id    = 0;
      name  = "Default Region";
      nodes = 0;
      faces = 0;
      delete[] m_node_id; m_node_id = nullptr;
      delete[] m_face_id; m_face_id = nullptr;
      m_domain = nullptr;
    }

    /// \brief Copy assignment operator
    gsCFX5Region& operator=(const gsCFX5Region& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Copy assigning gsCFX5Region\n";
#endif
          id    = other.id;
          name  = other.name;
          nodes = other.nodes;
          faces = other.nodes;

          delete[] m_node_id;
          m_node_id = new int[nodes];
          for (int i=0; i<nodes; i++)
            m_node_id[i] = other.m_node_id[i];

          delete[] m_face_id;
          m_face_id = new int[faces];
          for (int i=0; i<faces; i++)
            m_face_id[i] = other.m_face_id[i];

          m_domain = other.m_domain;
        }
      return *this;
    }

    /// \brief Move assignment operator
    gsCFX5Region& operator=(gsCFX5Region&& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Move assigning gsCFX5Region\n";
#endif
          id    = other.id;    other.id    = 0;
          name  = other.name;  other.name  = "Default Region";
          nodes = other.nodes; other.nodes = 0;
          faces = other.faces; other.faces = 0;

          delete[] m_node_id; m_node_id = other.m_node_id; other.m_node_id = nullptr;
          delete[] m_face_id; m_face_id = other.m_face_id; other.m_face_id = nullptr;

          m_domain = other.m_domain; other.m_domain  = nullptr;
        }
      return *this;
    }

    /// Returns swapped region
    gsCFX5Region& swap(gsCFX5Region& other)
    {
      std::swap(id,        other.id);
      std::swap(name,      other.name);
      std::swap(nodes,     other.nodes);
      std::swap(faces,     other.faces);
      std::swap(m_node_id, other.m_node_id);
      std::swap(m_face_id, other.m_face_id);
      std::swap(m_domain,  other.m_domain);
      return *this;
    }

  public:
    /// \brief Allocates memory for the list of node IDs and returns a
    /// pointer to it
    int* createNodeIDs(int nodes)
    {
      delete[] m_node_id;
      m_node_id = new int[nodes];

      return m_node_id;
    }

    /// \brief Allocates memory for the list of node IDs, updates the
    /// internal counter, and returns a pointer to the list of node
    /// IDs
    int* addNodeIDs(int nodes)
    {
      this->nodes = nodes;
      return createNodeIDs(nodes);
    }

    /// \brief Allocates memory for the list of face IDs and returns a
    /// pointer to it
    int* createFaceIDs(int facess)
    {
      delete[] m_face_id;
      m_face_id = new int[faces];

      return m_face_id;
    }

    /// \brief Allocates memory for the list of face IDs, updates the
    /// internal counter, and returns a pointer to the list of face
    /// IDs
    int* addFaceIDs(int faces)
    {
      this->faces = faces;
      return createFaceIDs(faces);
    }

    /// \brief Convert face ID from export to import format
    ///
    /// The CFX5 mesh import/export API uses different formats for
    /// storing face IDs. The internal storage of face IDs in the
    /// gsCFX5Region class adopts the cfxExport format. This methods
    /// converts face IDs stored in the cfxExport format to the
    /// cfxImport format.
    void convertFaceIDs();

  public:
    /// \brief Imports CFX5 region from CFX5 file
    void importFromCFX5(int id)
    {
      // Check that the region ID does not exceed the number of regions
      // present in the CFX5 file. Note that if the file has not been
      // initialized then the call to cfxExportRegionCount will return 0.
      if (id >= cfxExportRegionCount())
        GISMO_ERROR("Region ID exceeds number of regions in CFX5 file!");

      // Set the region ID
      this->id = id;

      // Import the region name from CFX5 file
      this->name = cfxExportRegionName(id+1);

      // Import the number of nodes from CFX5 file
      this->nodes = cfxExportRegionSize(id+1, cfxREG_NODES);

      // Import the list of node IDs from CFX5 file
      this->m_node_id = cfxCopy(cfxExportRegionList(id+1, cfxREG_NODES), this->nodes);

      // Import the number of faces from CFX5 file
      this->faces = cfxExportRegionSize(id+1, cfxREG_FACES);

      // Import the list of face IDs from CFX5 file
      this->m_face_id = cfxCopy(cfxExportRegionList(id+1, cfxREG_FACES), this->faces);

      // Free internal memory
      cfxExportRegionFree(id+1);
    }

    /// \brief Imports CFX5 region from ICEM ASCII file
    void importFromICEM(int id, std::ifstream& file)
    {
      if (!file.is_open())
        GISMO_ERROR("ICEM file is not ready for import!");

      // Set the region ID
      this->id = id;

      // Import number of face indices and region name from ICEM file
      file >> this->faces;
      file >> this->name;

      // Allocate memory for face IDs
      addFaceIDs(this->faces);

      // Import face IDs from ICEM file
      int element, face;
      for (int i=0; i<getFaceNumber(); i++)
        {
          file >> element;
          file >> face;
          this->m_face_id[i] = ((element << 3) | face);
        }

    }

    /// \brief Exports CFX5 region to CFX5 file
    void exportToCFX5()
    {
      // Check that the CFX5 file is ready
      if (cfxImportStatus() < 0)
        GISMO_ERROR("CFX5 file is not ready for import!");

      gsCFX5Region region = gsCFX5Region(*this);
      region.convertFaceIDs();

      // Export faces in region to CFX5 file
      if(cfxImportRegion(const_cast<char*>(region.getName().c_str()),
                         cfxImpREG_FACES, region.getFaceNumber(), (ID_t*)region.getFaceIDs())
         != (ID_t)region.getFaceNumber())
        GISMO_ERROR("An error occured while importing faces to CFX5 file!");
    }

    /// \brief Exports CFX5 region to ICEM ASCII file
    void exportToICEM(std::ofstream& file) const
    {
      if (!file.is_open())
        GISMO_ERROR("ICEM file is not ready for export!");

      // Export number of face IDs and region name to ICEM file
      file << this->faces << " " << this->name << std::endl;

      // Export face IDs to ICEM file
      for (auto faceidit = faceid_cbegin(); faceidit != faceid_cend(); faceidit++)
        file << cfxELEMNUM(*faceidit) << " " << cfxFACENUM(*faceidit) << std::endl;
    }

    /// \brief Returns the number of the region
    const int getRegion() const
    { return id; }

    /// \brief Sets the number of the region
    void setRegion(int id)
    { this->id = id; }

    /// \brief Returns the name of the region
    const std::string getName() const
    { return name; }

    /// \brief Sets the name of the region
    void setName(std::string name)
    { this->name = name; }

    /// \brief Returns the total number of nodes in the region
    const int getNodeNumber() const
    { return nodes; }

    /// \brief Returns the total number of faces in the region
    const int getFaceNumber() const
    { return faces; }

  public:
    /// \brief Returns constant pointer to the node ID array
    const int* getNodeIDs() const
    {
      return m_node_id;
    }

    /// \brief Returns pointer to node ID array
    int* getNodeIDs()
    {
      return m_node_id;
    }

    /// \brief Returns constant reference to the idx-th node ID of the element
    const int& getNodeID(const int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
      return m_node_id[idx];
    }

    /// \brief Returns reference to the idx-th node ID of the element
    int& getNodeID(const int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
      return m_node_id[idx];
    }

    /// \brief Returns constant reference to the idx-th node in the region
    const gsCFX5Node& getNode(int idx) const;

    /// \brief Returns reference to the idx-th node in the region
    gsCFX5Node& getNode(int idx);

  public:
    /// \brief Returns constant reference to the element associated with
    /// the idx-th face in the region
    const gsCFX5Element& getElement(int idx) const;

    /// \brief Returns reference to the element associated with the
    /// idx-th face in the region
    gsCFX5Element& getElement(int idx);

  public:
    /// \brief Returns constant pointer to the face ID array
    const int* getFaceIDs() const
    {
      return m_face_id;
    }

    /// \brief Returns pointer to face ID array
    int* getFaceIDs()
    {
      return m_face_id;
    }

    /// \brief Returns constant reference to the idx-th face ID of the region
    const int& getFaceID(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getFaceNumber()), "Invalid face index.");
      return m_face_id[idx];
    }

    /// \brief Returns reference to the idx-th face ID of the region
    int& getFaceID(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getFaceNumber()), "Invalid face index.");
      return m_face_id[idx];
    }

    /// \brief Returns idx-the face in the region
    gsCFX5Face getFace(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getFaceNumber()), "Invalid face index.");

      gsCFX5Element element = getElement(idx);
      switch (element.getElementType())
        {
        case gsCFX5ElementType::Tetrahedral :
          switch (cfxFACENUM(m_face_id[idx]))
            {
            case 1 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(1),
                                element.getNodeID(2));
            case 2 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(3),
                                element.getNodeID(1));
            case 3 :
              return gsCFX5Face(element.getNodeID(1),
                                element.getNodeID(3),
                                element.getNodeID(2));
            case 4 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(2),
                                element.getNodeID(3));
            default:
              GISMO_ERROR("Invalid face number!");
            }
          break;

        case gsCFX5ElementType::Pyramid :
          switch (cfxFACENUM(m_face_id[idx]))
            {
            case 1 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(3),
                                element.getNodeID(4));
            case 2 :
              return gsCFX5Face(element.getNodeID(1),
                                element.getNodeID(4),
                                element.getNodeID(2));
            case 3 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(4),
                                element.getNodeID(1));
            case 4 :
              return gsCFX5Face(element.getNodeID(2),
                                element.getNodeID(4),
                                element.getNodeID(3));
            case 5 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(1),
                                element.getNodeID(2),
                                element.getNodeID(3));
            default:
              GISMO_ERROR("Invalid face number!");
            }
          break;

        case gsCFX5ElementType::Prism :
          switch (cfxFACENUM(m_face_id[idx]))
            {
            case 1 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(2),
                                element.getNodeID(5),
                                element.getNodeID(3));
            case 2 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(3),
                                element.getNodeID(4),
                                element.getNodeID(1));
            case 3 :
              return gsCFX5Face(element.getNodeID(1),
                                element.getNodeID(4),
                                element.getNodeID(5),
                                element.getNodeID(2));
            case 4 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(1),
                                element.getNodeID(2));
            case 5 :
              return gsCFX5Face(element.getNodeID(3),
                                element.getNodeID(5),
                                element.getNodeID(4));
            default:
              GISMO_ERROR("Invalid face number!");
            }
          break;


        case gsCFX5ElementType::Hexahedral :
          switch (cfxFACENUM(m_face_id[idx]))
            {
            case 1 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(2),
                                element.getNodeID(6),
                                element.getNodeID(4));
            case 2 :
              return gsCFX5Face(element.getNodeID(1),
                                element.getNodeID(5),
                                element.getNodeID(7),
                                element.getNodeID(3));
            case 3 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(4),
                                element.getNodeID(5),
                                element.getNodeID(1));
            case 4 :
              return gsCFX5Face(element.getNodeID(2),
                                element.getNodeID(3),
                                element.getNodeID(7),
                                element.getNodeID(6));
            case 5 :
              return gsCFX5Face(element.getNodeID(0),
                                element.getNodeID(1),
                                element.getNodeID(3),
                                element.getNodeID(2));
            case 6 :
              return gsCFX5Face(element.getNodeID(4),
                                element.getNodeID(6),
                                element.getNodeID(7),
                                element.getNodeID(5));
            default:
              GISMO_ERROR("Invalid face number!");
            }
          break;

        default:
          GISMO_ERROR("Invalid element type!");
        }
    }

    /// \brief Returns hash of the region
    size_t getHash() const
    {
      std::string region = util::to_string(getRegion())+getName()+
        util::to_string(getNodeNumber())+util::to_string(getFaceNumber());

      size_t hash = std::hash<std::string>{}(region);

      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(nodeit->getHash()));

      for (auto faceit = face_cbegin(); faceit != face_cend(); faceit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*faceit).getHash()));

      return hash;
    }

    /// \brief Prints the region
    std::ostream& print(std::ostream& os) const
    {
      os << "CFX5Region (" << getRegion() << "): " << getName() << "\n"
         << std::setw(20)  << getNodeNumber()      << "    nodes\n"
         << std::setw(20)  << getFaceNumber()      << "    faces";
      return os;
    }

    /// \brief Prints the region as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);

      os << "\n========== NODES IN REGION ==========\n";
      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        os << (*nodeit).detail() << std::endl;

      os << "\n========== FACES IN REGION ==========\n";
      for (auto faceit = face_cbegin(); faceit != face_cend(); faceit++)
        os << (*faceit).detail() << std::endl;

      return os.str();
    }

    /// \brief Prints the difference between the region and the other
    /// region as a string
    std::string diff(const gsCFX5Region& other) const
    {
      std::ostringstream os;

      if (getRegion() != other.getRegion())
        os << "CFX5Region region\n< "
           << getRegion()
           << "\n---\n"
           << "> "
           << other.getRegion()
           << "\n";

      if (getNodeNumber() != other.getNodeNumber())
        os << "CFX5Region nodes\n< "
           << getNodeNumber()
           << "\n---\n"
           << "> "
           << other.getNodeNumber()
           << "\n";

      if (getFaceNumber() != other.getFaceNumber())
        os << "CFX5Region faces\n< "
           << getFaceNumber()
           << "\n---\n"
           << "> "
           << other.getFaceNumber()
           << "\n";

      os << "\n========== NODES IN REGION ==========\n";
      for (auto nodeit = node_cbegin(), other_nodeit = other.node_cbegin();
           nodeit != node_cend() && other_nodeit != other.node_cend();
           nodeit++, other_nodeit++)
        os << (*nodeit).diff(*other_nodeit);

      os << "\n========== FACES IN REGION ==========\n";
      for (auto faceit = face_cbegin(), other_faceit = other.face_cbegin();
           faceit != face_cend() && other_faceit != other.face_cend();
           faceit++, other_faceit++)
        os << (*faceit).diff(*other_faceit);

      return os.str();
    }

  public:
    typedef gsCFX5GenericIterator<int>             nodeid_iterator;
    typedef gsCFX5GenericConstIterator<int>        nodeid_const_iterator;
    typedef gsCFX5GenericReverseIterator<int>      nodeid_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<int> nodeid_const_reverse_iterator;

    /// Returns an iterator to the beginning of the node index list
    nodeid_iterator nodeid_begin()
    {
      return nodeid_iterator(m_node_id);
    }

    /// Returns a const-iterator to the beginning of the node index list
    nodeid_const_iterator nodeid_cbegin() const
    {
      return nodeid_const_iterator(m_node_id);
    }

    /// Returns a reverse-iterator to the beginning of the node index list
    nodeid_reverse_iterator nodeid_rbegin()
    {
      return nodeid_reverse_iterator(m_node_id + getNodeNumber());
    }

    /// Returns a const-reverse-iterator to the beginning of the node index list
    nodeid_const_reverse_iterator nodeid_crbegin() const
    {
      return nodeid_const_reverse_iterator(m_node_id + getNodeNumber());
    }

    /// Returns an iterator to the end of the node index list
    nodeid_iterator nodeid_end()
    {
      return nodeid_iterator(m_node_id + getNodeNumber());
    }

    /// Returns a const-iterator to the end of the node index list
    nodeid_const_iterator nodeid_cend() const
    {
      return nodeid_const_iterator(m_node_id + getNodeNumber());
    }

    /// Returns a reverse-iterator to the end of the node index list
    nodeid_reverse_iterator nodeid_rend()
    {
      return nodeid_reverse_iterator(m_node_id);
    }

    /// Returns a const-reverse-iterator to the end of the node index list
    nodeid_const_reverse_iterator nodeid_crend() const
    {
      return nodeid_const_reverse_iterator(m_node_id);
    }

  public:
    typedef gsCFX5NodeInRegionIterator             node_iterator;
    typedef gsCFX5NodeInRegionConstIterator        node_const_iterator;
    typedef gsCFX5NodeInRegionReverseIterator      node_reverse_iterator;
    typedef gsCFX5NodeInRegionConstReverseIterator node_const_reverse_iterator;

    /// Returns an iterator to the beginning of the node list
    node_iterator node_begin();

    /// Returns a const-iterator to the beginning of the node list
    node_const_iterator node_cbegin() const;

    /// Returns a reverse-iterator to the beginning of the node list
    node_reverse_iterator node_rbegin();

    /// Returns a const-reverse-iterator to the beginning of the node list
    node_const_reverse_iterator node_crbegin() const;

    /// Returns an iterator to the end of the node list
    node_iterator node_end();

    /// Returns a const-iterator to the end of the node list
    node_const_iterator node_cend() const;

    /// Returns a reverse-iterator to the end of the node list
    node_reverse_iterator node_rend();

    /// Returns a const-reverse-iterator to the end of the node list
    node_const_reverse_iterator node_crend() const;

  public:
    typedef gsCFX5GenericIterator<int>             faceid_iterator;
    typedef gsCFX5GenericConstIterator<int>        faceid_const_iterator;
    typedef gsCFX5GenericReverseIterator<int>      faceid_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<int> faceid_const_reverse_iterator;

    /// Returns an iterator to the beginning of the face index list
    faceid_iterator faceid_begin()
    {
      return faceid_iterator(m_face_id);
    }

    /// Returns a const-iterator to the beginning of the face index list
    faceid_const_iterator faceid_cbegin() const
    {
      return faceid_const_iterator(m_face_id);
    }

    /// Returns a reverse-iterator to the beginning of the face index list
    faceid_reverse_iterator faceid_rbegin()
    {
      return faceid_reverse_iterator(m_face_id + getFaceNumber());
    }

    /// Returns a const-reverse-iterator to the beginning of the face index list
    faceid_const_reverse_iterator faceid_crbegin() const
    {
      return faceid_const_reverse_iterator(m_face_id + getFaceNumber());
    }

    /// Returns an iterator to the end of the face index list
    faceid_iterator faceid_end()
    {
      return faceid_iterator(m_face_id + getFaceNumber());
    }

    /// Returns a const-iterator to the end of the face index list
    faceid_const_iterator faceid_cend() const
    {
      return faceid_const_iterator(m_face_id + getFaceNumber());
    }

    /// Returns a reverse-iterator to the end of the face index list
    faceid_reverse_iterator faceid_rend()
    {
      return faceid_reverse_iterator(m_face_id);
    }

    /// Returns a const-reverse-iterator to the end of the face index list
    faceid_const_reverse_iterator faceid_crend() const
    {
      return faceid_const_reverse_iterator(m_face_id);
    }

  public:
    typedef gsCFX5FaceInRegionConstIterator        face_const_iterator;
    typedef gsCFX5FaceInRegionConstReverseIterator face_const_reverse_iterator;

    /// Returns a const-iterator to the beginning of the face list
    face_const_iterator face_cbegin() const
    {
      return face_const_iterator(0, this);
    }

    /// Returns a const-reverse-iterator to the beginning of the face list
    face_const_reverse_iterator face_crbegin() const
    {
      return face_const_reverse_iterator(getFaceNumber(), this);
    }

    /// Returns a const-iterator to the end of the face list
    face_const_iterator face_cend() const
    {
      return face_const_iterator(getFaceNumber(), this);
    }

    /// Returns a const-reverse-iterator to the end of the face list
    face_const_reverse_iterator face_crend() const
    {
      return face_const_reverse_iterator(0, this);
    }

  protected:
    // Region ID
    int id;

    // Region name
    std::string name;

    // Number of nodes in the region
    int nodes;

    // Number of faces in the region
    int faces;

    // List of node IDs in the region
    int* m_node_id = nullptr;

    // List of face IDs in the region
    int* m_face_id = nullptr;

    // Pointer to the underlying CFX5 domain
    gsCFX5Domain* m_domain = nullptr;
  };

  /// \brief Output CFX5 region object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Region& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 volume.
  */
  class gsCFX5Volume
  {
  public:
    /// \brief Constructor (default)
    gsCFX5Volume()
      :
      id{0},
      name{"Default Volume"},
      nodes{0},
      elements{0},
      m_node_id{nullptr},
      m_element_id{nullptr},
      m_domain{nullptr}
    {}

    /// \brief Constructor (with pointer to underlying domain)
    gsCFX5Volume(gsCFX5Domain* m_domain)
      :
      id{0},
      name{"Default Volume"},
      nodes{0},
      elements{0},
      m_node_id{nullptr},
      m_element_id{nullptr},
      m_domain{m_domain}
    {}

    /// \brief Constructor (copy)
    gsCFX5Volume(const gsCFX5Volume& other)
      :
      id{other.id},
      name{other.name},
      nodes{other.nodes},
      elements{other.elements},
      m_node_id{new int[other.nodes]},
      m_element_id{new int[other.elements]},
      m_domain{other.m_domain}
    {
#if !defined(NDEBUG)
      gsDebug << "Copying gsCFX5Volume\n";
#endif
      for (int i=0; i<nodes; i++)
        m_node_id[i] = other.m_node_id[i];

      for (int i=0; i<elements; i++)
        m_element_id[i] = other.m_element_id[i];
    }

    /// \brief Constructor (move)
    gsCFX5Volume(gsCFX5Volume&& other)
      :
      id{other.id},
      name{other.name},
      nodes{other.nodes},
      elements{other.elements},
      m_node_id{other.m_node_id},
      m_element_id{other.m_element_id},
      m_domain{other.m_domain}
    {
#if !defined(NDEBUG)
      gsDebug << "Moving gsCFX5Volume\n";
#endif
      other.id           = 0;
      other.name         = "Default Volume";
      other.nodes        = 0;
      other.elements     = 0;
      other.m_node_id    = nullptr;
      other.m_element_id = nullptr;
      other.m_domain     = nullptr;
    }

    /// \brief Destructor
    ~gsCFX5Volume()
    {
      id       = 0;
      name     = "Default Volume";
      nodes    = 0;
      elements = 0;
      delete[] m_node_id;    m_node_id = nullptr;
      delete[] m_element_id; m_element_id = nullptr;
      m_domain = nullptr;
    }

    /// \brief Copy assignment operator
    gsCFX5Volume& operator=(const gsCFX5Volume& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Copy assigning gsCFX5Volume\n";
#endif
          id       = other.id;
          name     = other.name;
          nodes    = other.nodes;
          elements = other.elements;

          delete[] m_node_id;
          m_node_id = new int[nodes];
          for (int i=0; i<nodes; i++)
            m_node_id[i] = other.m_node_id[i];

          delete[] m_element_id;
          m_element_id = new int[elements];
          for (int i=0; i<elements; i++)
            m_element_id[i] = other.m_element_id[i];

          m_domain = other.m_domain;
        }
      return *this;
    }

    /// \brief Move assignment operator
    gsCFX5Volume& operator=(gsCFX5Volume&& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Move assigning gsCFX5Volume\n";
#endif
          id       = other.id;       other.id       = 0;
          name     = other.name;     other.name     = "Default Volume";
          nodes    = other.nodes;    other.nodes    = 0;
          elements = other.elements; other.elements = 0;

          delete[] m_node_id;
          m_node_id = other.m_node_id; other.m_node_id = nullptr;

          delete[] m_element_id;
          m_element_id = other.m_element_id; other.m_element_id = nullptr;

          m_domain = other.m_domain; other.m_domain  = nullptr;
        }
      return *this;
    }

    /// Returns swapped volume
    gsCFX5Volume& swap(gsCFX5Volume& other)
    {
      std::swap(id,           other.id);
      std::swap(name,         other.name);
      std::swap(nodes,        other.nodes);
      std::swap(elements,     other.elements);
      std::swap(m_node_id,    other.m_node_id);
      std::swap(m_element_id, other.m_element_id);
      std::swap(m_domain,     other.m_domain);
      return *this;
    }

  public:
    /// \brief Allocates memory for the list of node IDs and returns a
    /// pointer to it
    int* createNodeIDs(int nodes)
    {
      delete[] m_node_id;
      m_node_id = new int[nodes];

      return m_node_id;
    }

    /// \brief Allocates memory for the list of node IDs, updates the
    /// internal counter, and returns a pointer to the list of node
    /// IDs
    int* addNodeIDs(int nodes)
    {
      this->nodes = nodes;
      return createNodeIDs(nodes);
    }

    /// \brief Allocates memory for the list of element IDs and
    /// returns a pointer to it
    int* createElementIDs(int elements)
    {
      delete[] m_element_id;
      m_element_id = new int[elements];

      return m_element_id;
    }

    /// \brief Allocates memory for the list of element IDs, updates
    /// the internal counter, and returns a pointer to the list of
    /// element IDs
    int* addElementIDs(int elements)
    {
      this->elements = elements;
      return createElementIDs(elements);
    }

  public:
    /// \brief Imports CFX5 volume from CFX5 file
    void importFromCFX5(int id)
    {
      // Check that the volume ID does not exceed the number of volumes
      // present in the CFX5 file. Note that if the file has not been
      // initialized then the call to cfxExportVolumeCount will return 0.
      if (id >= cfxExportVolumeCount())
        GISMO_ERROR("Volume ID exceeds number of volumes in CFX5 file!");

      // Set the volume ID
      this->id = id;

      // Import the volume name from CFX5 file
      this->name = cfxExportVolumeName(id+1);

      // Import the number of nodes from CFX5 file
      this->nodes = cfxExportVolumeSize(id+1, cfxVOL_NODES);

      // Import the list of node IDs from CFX5 file
      this->m_node_id = cfxCopy(cfxExportVolumeList(id+1, cfxVOL_NODES), this->nodes);

      // Import the number of elements from CFX5 file
      this->elements = cfxExportVolumeSize(id+1, cfxVOL_ELEMS);

      // Import the list of element IDs from CFX5 file
      this->m_element_id = cfxCopy(cfxExportVolumeList(id+1, cfxVOL_ELEMS), this->elements);

      // Free internal memory
      cfxExportVolumeFree(id+1);
    }

    /// \brief Imports CFX5 volume from ICEM ASCII file
    void importFromICEM(int id, std::ifstream& file)
    {
      if (!file.is_open())
        GISMO_ERROR("ICEM file is not ready for import!");

      // Set the volume ID
      this->id = id;

      // Import number of element indices and volume name
      int elements;
      file >> elements;
      file >> this->name;

      // Allocate memory for element IDs
      addElementIDs(elements);

      // Import element IDs
      for (int i=0; i<getElementNumber(); i++)
        file >> m_element_id[i];
    }

    /// \brief Exports CFX5 volume to CFX5 file
    void exportToCFX5() const
    {
      // Check that the CFX5 file is ready
      if (cfxImportStatus() < 0)
        GISMO_ERROR("CFX5 file is not ready for import!");

      // Export elements to CFX5 file
      if (cfxImportRegion(const_cast<char*>(getName().c_str()),
                          cfxImpREG_ELEMS, getElementNumber(), (ID_t*)getElementIDs())
          != (ID_t)getElementNumber())
        GISMO_ERROR("An error occured while importing volumes to CFX5 file!");
    }

    /// \brief Exports CFX5 volume to ICEM ASCII file
    void exportToICEM(std::ofstream& file) const
    {
      if (!file.is_open())
        GISMO_ERROR("ICEM file is not ready for export!");

      // Export number of element IDs and volume name to ICEM file
      file << this->elements << " " << this->name << std::endl;

      // Export element IDs to ICEM file
      for (auto elemidit = elementid_cbegin(); elemidit != elementid_cend(); elemidit++)
        file << *elemidit << std::endl;
    }

    /// \brief Returns the number of the volume
    const int getVolume() const
    { return id; }

    /// \brief Sets the number of the volume
    void setVolume(int id)
    { this->id = id; }

    /// \brief Returns the name of the volume
    const std::string getName() const
    { return name; }

    /// \brief Sets the name of the volume
    void setName(const std::string& name)
    { this->name = name; }

    /// \brief Returns the total number of nodes in the volume
    const int getNodeNumber() const
    { return nodes; }

    /// \brief Returns the total number of elements in the volume
    const int getElementNumber() const
    { return elements; }

  public:
    /// \brief Returns constant pointer to the node ID array
    const int* getNodeIDs() const
    {
      return m_node_id;
    }

    /// \brief Returns pointer to the node ID array
    int* getNodeIDs()
    {
      return m_node_id;
    }

    /// \brief Returns constant reference to the idx-th node
    const gsCFX5Node& getNode(int idx) const;

    /// \brief Returns reference to the idx-th node
    gsCFX5Node& getNode(int idx);

    /// \brief Returns hash of the volume
    size_t getHash() const
    {
      std::string volume = util::to_string(getVolume())+getName()+
        util::to_string(getNodeNumber())+util::to_string(getElementNumber());

      size_t hash = std::hash<std::string>{}(volume);

      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(nodeit->getHash()));

      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(elemit->getHash()));

      return hash;
      }

    /// \brief Prints the volume
    std::ostream& print(std::ostream& os) const
    {
      os << "CFX5Volume (" << getVolume() << "): " << getName() << "\n"
         << std::setw(20)  << getNodeNumber()      << "    nodes\n"
         << std::setw(20)  << getElementNumber()   << "    elements";
      return os;
    }

    /// \brief Prints the volume as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);

      os << "\n========== NODES IN VOLUME ==========\n";
      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        os << (*nodeit).detail() << std::endl;

      os << "\n========== ELEMENTS IN VOLUME ==========\n";
      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++)
        os << (*elemit).detail() << std::endl;

      return os.str();
    }

    /// \brief Prints the difference between the volume and the other
    /// volume as a string
    std::string diff(const gsCFX5Volume& other) const
    {
      std::ostringstream os;

      if (getVolume() != other.getVolume())
        os << "CFX5Volume volume\n< "
           << getVolume()
           << "\n---\n"
           << "> "
           << other.getVolume()
           << "\n";

      if (getNodeNumber() != other.getNodeNumber())
        os << "CFX5Volume nodes\n< "
           << getNodeNumber()
           << "\n---\n"
           << "> "
           << other.getNodeNumber()
           << "\n";

      if (getElementNumber() != other.getElementNumber())
        os << "CFX5Volume elements\n< "
           << getElementNumber()
           << "\n---\n"
           << "> "
           << other.getElementNumber()
           << "\n";

      os << "\n========== NODES IN REGION ==========\n";
      for (auto nodeit = node_cbegin(), other_nodeit = other.node_cbegin();
           nodeit != node_cend() && other_nodeit != other.node_cend();
           nodeit++, other_nodeit++)
        os << (*nodeit).diff(*other_nodeit);

      os << "\n========== ELEMENTS IN REGION ==========\n";
      for (auto elemit = element_cbegin(), other_elemit = other.element_cbegin();
           elemit != element_cend() && other_elemit != other.element_cend();
           elemit++, other_elemit++)
        os << (*elemit).diff(*other_elemit);

      return os.str();
    }

  public:
    typedef gsCFX5GenericIterator<int>             nodeid_iterator;
    typedef gsCFX5GenericConstIterator<int>        nodeid_const_iterator;
    typedef gsCFX5GenericReverseIterator<int>      nodeid_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<int> nodeid_const_reverse_iterator;

    /// Returns an iterator to the beginning of the node index list
    nodeid_iterator nodeid_begin()
    {
      return nodeid_iterator(m_node_id);
    }

    /// Returns a const-iterator to the beginning of the node index list
    nodeid_const_iterator nodeid_cbegin() const
    {
      return nodeid_const_iterator(m_node_id);
    }

    /// Returns a reverse-iterator to the beginning of the node index list
    nodeid_reverse_iterator nodeid_rbegin()
    {
      return nodeid_reverse_iterator(m_node_id + getNodeNumber());
    }

    /// Returns a const-reverse-iterator to the beginning of the node index list
    nodeid_const_reverse_iterator nodeid_crbegin() const
    {
      return nodeid_const_reverse_iterator(m_node_id + getNodeNumber());
    }

    /// Returns an iterator to the end of the node index list
    nodeid_iterator nodeid_end()
    {
      return nodeid_iterator(m_node_id + getNodeNumber());
    }

    /// Returns a const-iterator to the end of the node index list
    nodeid_const_iterator nodeid_cend() const
    {
      return nodeid_const_iterator(m_node_id + getNodeNumber());
    }

    /// Returns a reverse-iterator to the end of the node index list
    nodeid_reverse_iterator nodeid_rend()
    {
      return nodeid_reverse_iterator(m_node_id);
    }

    /// Returns a const-reverse-iterator to the end of the node index list
    nodeid_const_reverse_iterator nodeid_crend() const
    {
      return nodeid_const_reverse_iterator(m_node_id);
    }

  public:
    typedef gsCFX5NodeInVolumeIterator             node_iterator;
    typedef gsCFX5NodeInVolumeConstIterator        node_const_iterator;
    typedef gsCFX5NodeInVolumeReverseIterator      node_reverse_iterator;
    typedef gsCFX5NodeInVolumeConstReverseIterator node_const_reverse_iterator;

    /// Returns an iterator to the beginning of the node list
    node_iterator node_begin();

    /// Returns a const-iterator to the beginning of the node list
    node_const_iterator node_cbegin() const;

    /// Returns a reverse-iterator to the beginning of the node list
    node_reverse_iterator node_rbegin();

    /// Returns a const-reverse-iterator to the beginning of the node list
    node_const_reverse_iterator node_crbegin() const;

    /// Returns an iterator to the end of the node list
    node_iterator node_end();

    /// Returns a const-iterator to the end of the node list
    node_const_iterator node_cend() const;

    /// Returns a reverse-iterator to the end of the node list
    node_reverse_iterator node_rend();

    /// Returns a const-reverse-iterator to the end of the node list
    node_const_reverse_iterator node_crend() const;

  public:
    /// \brief Returns constant pointer to the element ID array
    const int* getElementIDs() const
    {
      return m_element_id;
    }

    /// \brief Returns pointer to element ID array
    int* getElementIDs()
    {
      return m_element_id;
    }

  public:
    typedef gsCFX5GenericIterator<int>             elementid_iterator;
    typedef gsCFX5GenericConstIterator<int>        elementid_const_iterator;
    typedef gsCFX5GenericReverseIterator<int>      elementid_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<int> elementid_const_reverse_iterator;

    /// Returns an iterator to the beginning of the element index list
    elementid_iterator elementid_begin()
    {
      return elementid_iterator(m_element_id);
    }

    /// Returns a const-iterator to the beginning of the element index list
    elementid_const_iterator elementid_cbegin() const
    {
      return elementid_const_iterator(m_element_id);
    }

    /// Returns a reverse-iterator to the beginning of the element index list
    elementid_reverse_iterator elementid_rbegin()
    {
      return elementid_reverse_iterator(m_element_id + getElementNumber());
    }

    /// Returns a const-reverse-iterator to the beginning of the element index list
    elementid_const_reverse_iterator elementid_crbegin() const
    {
      return elementid_const_reverse_iterator(m_element_id + getElementNumber());
    }

    /// Returns an iterator to the end of the element index list
    elementid_iterator elementid_end()
    {
      return elementid_iterator(m_element_id + getElementNumber());
    }

    /// Returns a const-iterator to the end of the element index list
    elementid_const_iterator elementid_cend() const
    {
      return elementid_const_iterator(m_element_id + getElementNumber());
    }

    /// Returns a reverse-iterator to the end of the element index list
    elementid_reverse_iterator elementid_rend()
    {
      return elementid_reverse_iterator(m_element_id);
    }

    /// Returns a const-reverse-iterator to the end of the element index list
    elementid_const_reverse_iterator elementid_crend() const
    {
      return elementid_const_reverse_iterator(m_element_id);
    }

  public:
    typedef gsCFX5ElementInVolumeIterator             element_iterator;
    typedef gsCFX5ElementInVolumeConstIterator        element_const_iterator;
    typedef gsCFX5ElementInVolumeReverseIterator      element_reverse_iterator;
    typedef gsCFX5ElementInVolumeConstReverseIterator element_const_reverse_iterator;

    /// Returns an iterator to the beginning of the element list
    element_iterator element_begin();

    /// Returns a const-iterator to the beginning of the element list
    element_const_iterator element_cbegin() const;

    /// Returns a reverse-iterator to the beginning of the element list
    element_reverse_iterator element_rbegin();

    /// Returns a const-reverse-iterator to the beginning of the element list
    element_const_reverse_iterator element_crbegin() const;

    /// Returns an iterator to the end of the element list
    element_iterator element_end();

    /// Returns a const-iterator to the end of the element list
    element_const_iterator element_cend() const;

    /// Returns a reverse-iterator to the end of the element list
    element_reverse_iterator element_rend();

    /// Returns a const-reverse-iterator to the end of the element list
    element_const_reverse_iterator element_crend() const;

  private:
    // Volume ID
    int id;

    // Volume name
    std::string name;

    // Number of nodes in the volume
    int nodes;

    // Number of elements in the volume
    int elements;

    // List of node IDs in the volume
    int* m_node_id = nullptr;

    // List of element IDs in the volume
    int* m_element_id = nullptr;

    // Pointer to the underlying CFX5 domain
    gsCFX5Domain* m_domain = nullptr;
  };

  /// \brief Output CFX5 volume object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Volume& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 variable.
  */
  class gsCFX5Variable
  {
  public:
    /// \brief Constructor (default)
    gsCFX5Variable()
      :
      id{0},
      alias{"Default Variable"},
      name{"Default Variable"},
      length{0},
      dimension{0},
      bdrflag{0},
      m_value{nullptr}
    {
#if !defined(NDEBUG)
      gsDebug << "Constructing gsCFX5Variable\n";
#endif
    }

    /// \brief Constructor (copy)
    gsCFX5Variable(const gsCFX5Variable& other)
      :
      id{other.id},
                 alias{other.alias},
                 name{other.name},
                 length{other.length},
                 dimension{other.dimension},
                 bdrflag{other.bdrflag},
                 m_value{new float[other.length]}
    {
#if !defined(NDEBUG)
      gsDebug << "Copying gsCFX5Variable\n";
#endif
      for (int i=0; i<length; i++)
        m_value[i] = other.m_value[i];
    }

    /// \brief Constructor (move)
    gsCFX5Variable(gsCFX5Variable&& other)
      :
      id{other.id},
                 alias{other.alias},
                 name{other.name},
                 length{other.length},
                 dimension{other.dimension},
                 bdrflag{other.bdrflag},
                 m_value{other.m_value}
    {
#if !defined(NDEBUG)
      gsDebug << "Moving gsCFX5Variable\n";
#endif
      other.id        = 0;
      other.alias     = "Default Variable";
      other.name      = "Default Variable";
      other.length    = 0;
      other.dimension = 0;
      other.bdrflag   = 0;
      other.m_value   = nullptr;
    }

    /// \brief Destructor
    ~gsCFX5Variable()
    {
      id        = 0;
      alias     = "Default Variable";
      name      = "Default Variable";
      length    = 0;
      dimension = 0;
      bdrflag   = 0;
      delete[] m_value; m_value = nullptr;
    }

    /// \brief Copy assignment operator
    gsCFX5Variable& operator=(const gsCFX5Variable& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Copy assigning gsCFX5Variable\n";
#endif
          id        = other.id;
          alias     = other.alias;
          name      = other.name;
          length    = other.length;
          dimension = other.dimension;
          bdrflag   = other.bdrflag;

          delete[] m_value;
          m_value = new float[length];
          for (int i=0; i<length; i++)
            m_value[i] = other.m_value[i];
        }
      return *this;
    }

    /// \brief Move assignment operator
    gsCFX5Variable& operator=(gsCFX5Variable&& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Move assigning gsCFX5Variable\n";
#endif
          id        = other.id;        other.id        = 0;
          alias     = other.alias;     other.alias     = "Default Variable";
          name      = other.name;      other.name      = "Default Variable";
          length    = other.length;    other.length    = 0;
          dimension = other.dimension; other.dimension = 0;
          bdrflag   = other.bdrflag;   other.bdrflag   = 0;

          delete[] m_value;
          m_value = other.m_value; other.m_value = nullptr;
        }
      return *this;
    }

    /// \brief Returns swapped variable
    gsCFX5Variable& swap(gsCFX5Variable& other)
    {
      std::swap(id,        other.id);
      std::swap(alias,     other.alias);
      std::swap(name,      other.name);
      std::swap(length,    other.length);
      std::swap(dimension, other.dimension);
      std::swap(bdrflag,   other.bdrflag);
      std::swap(m_value,   other.m_value);
      return *this;
    }

  public:
    typedef gsCFX5GenericIterator<float>             value_iterator;
    typedef gsCFX5GenericConstIterator<float>        value_const_iterator;
    typedef gsCFX5GenericReverseIterator<float>      value_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<float> value_const_reverse_iterator;

    /// Returns an iterator to the beginning of the value list
    value_iterator begin()
    {
      return value_iterator(m_value);
    }

    /// Returns a const-iterator to the beginning of the value list
    value_const_iterator cbegin() const
    {
      return value_const_iterator(m_value);
    }

    /// Returns a reverse-iterator to the beginning of the value list
    value_reverse_iterator rbegin()
    {
      return value_reverse_iterator(m_value + getLength());
    }

    /// Returns a const-reverse-iterator to the beginning of the value list
    value_const_reverse_iterator crbegin() const
    {
      return value_const_reverse_iterator(m_value + getLength());
    }

    /// Returns an iterator to the end of the value list
    value_iterator end()
    {
      return value_iterator(m_value + getLength());
    }

    /// Returns a const-iterator to the end of the value list
    value_const_iterator cend() const
    {
      return value_const_iterator(m_value + getLength());
    }

    /// Returns a reverse-iterator to the end of the value list
    value_reverse_iterator rend()
    {
      return value_reverse_iterator(m_value);
    }

    /// Returns a const-reverse-iterator to the end of the value list
    value_const_reverse_iterator crend() const
    {
      return value_const_reverse_iterator(m_value);
    }

  public:
    /// \brief Imports CFX5 variable from CFX5 file
    void importFromCFX5(int id)
    {
      // Check that the variable ID does not exceed the number of variables
      // present in the CFX5 file. Note that if the file has not been
      // initialized then the call to cfxExportVariableCount will return 0.
      if (id >= cfxExportVariableCount(0))
        GISMO_ERROR("Variable ID exceeds number of variables in CFX5 file!");

      // Set the variable ID
      this->id = id;

      // Import the variable alias from CFX5 file
      this->alias = cfxExportVariableName(id+1, 0);

      // Import the variable name from CFX5 file
      this->name = cfxExportVariableName(id+1, 1);

      // Import the dimension and the length of the variable the CFX5 file
      int status = cfxExportVariableSize(id+1, &this->dimension, &this->length, &this->bdrflag);
      GISMO_ASSERT(status != 0, "An error occured while importing variable from CFX5 file.");
      GISMO_UNUSED(status);

      // Import the variable data from CFX5 file
      this->m_value = cfxCopy(cfxExportVariableList(id+1, 1), this->length*this->dimension);

      // Free internal memory
      cfxExportVariableFree(id+1);
    }

    /// \brief Returns the number of the variable
    const int getVariable() const
    { return id; }

    /// \brief Sets the number of the variable
    void setVariable(int id)
    { this->id = id; }

    /// \brief Returns the alias of the variable
    const std::string& getAlias() const
    { return alias; }

    /// \brief Sets the alias of the variable
    void setAlias(const std::string& name)
    { this->alias = alias; }

    /// \brief Returns the name of the variable
    const std::string& getName() const
    { return name; }

    /// \brief Sets the name of the variable
    void setName(const std::string& name)
    { this->name = name; }

    /// \brief Returns the dimension of the variable
    const int getDimension() const
    { return dimension; }

    /// \brief Sets the dimension of the variable
    void setDimension(int dimension)
    { this->dimension = dimension; }

    /// \brief Returns the length of the variable
    const int getLength() const
    { return length; }

    /// \brief Returns the boundary flag of the variable
    const int getBdrFlag() const
    { return bdrflag; }

    /// \brief Sets the boundary flag of the variable
    void setBdrFlag(bool bdrflag)
    { this->bdrflag = bdrflag; }

  public:
    /// \brief Returns constant pointer to the value array
    const float* getValues() const
    {
      return m_value;
    }

    /// \brief Returns pointer to the value array
    float* getValues()
    {
      return m_value;
    }

    /// \brief Returns constant reference to the idx-th value of the variable
    const float& getValue(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getLength()), "Invalid value index.");
      return m_value[idx];
    }

    /// \brief Returns reference to the idx-th value of the variable
    float& getValue(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getLength()), "Invalid value index.");
      return m_value[idx];
    }

    /// \brief Returns hash of the variable
    size_t getHash() const
    {
      std::string variable = util::to_string(getVariable())+getName()+getAlias()+
        util::to_string(getLength())+util::to_string(getDimension())+
        util::to_string(getBdrFlag());

      size_t hash = std::hash<std::string>{}(variable);

      for (auto valueit = cbegin(); valueit != cend(); valueit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(*valueit));

      return hash;
    }

    /// \brief Prints the variable
    std::ostream& print(std::ostream& os) const
    {
      os << "CFX5Variable (" << getVariable()   << "): "
         << getName() << " [" << getAlias()     << "]\n"
         << std::setw(20)     << getLength()    << "    length\n"
         << std::setw(20)     << getDimension() << "    dimension\n"
         << std::setw(20)     << getBdrFlag()   << "    bdrflag";
      return os;
    }

    /// \brief Prints the variable as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);

      os << "\n========== VALUE IN VARIABLE ==========\n";
      for (auto valueit = cbegin(); valueit != cend(); valueit++)
        os << *valueit << std::endl;

      return os.str();
    }

    /// \brief Prints the difference between the variable and the other
    /// variable as a string
    std::string diff(const gsCFX5Variable& other) const
    {
      std::ostringstream os;

      if (getVariable() != other.getVariable())
        os << "CFX5Variable variable\n< "
           << getVariable()
           << "\n---\n"
           << "> "
           << other.getVariable()
           << "\n";

      if (getName() != other.getName())
        os << "CFX5Variable name\n< "
           << getName()
           << "\n---\n"
           << "> "
           << other.getName()
           << "\n";

      if (getAlias() != other.getAlias())
        os << "CFX5Variable alias\n< "
           << getAlias()
           << "\n---\n"
           << "> "
           << other.getAlias()
           << "\n";

      if (getLength() != other.getLength())
        os << "CFX5Variable length\n< "
           << getLength()
           << "\n---\n"
           << "> "
           << other.getLength()
           << "\n";

      if (getDimension() != other.getDimension())
        os << "CFX5Variable dimension\n< "
           << getDimension()
           << "\n---\n"
           << "> "
           << other.getDimension()
           << "\n";

      if (getBdrFlag() != other.getBdrFlag())
        os << "CFX5Variable bdrflag\n< "
           << getBdrFlag()
           << "\n---\n"
           << "> "
           << other.getBdrFlag()
           << "\n";

      os << "\n========== VALUE IN VARIABLE ==========\n";
      for (auto valit = cbegin(), other_valit = other.cbegin();
           valit != cend() && other_valit != other.cend();
           valit++, other_valit++)
        if (*valit != *other_valit)
          os << "CFX5Variable value\n< "
             << *valit << "\n---\n> " << *other_valit << "\n";

      return os.str();
    }

  private:
    // Variable ID
    int id;

    // Variable alias
    std::string alias;

    // Variable name
    std::string name;

    // Length of the variable
    int length;

    // Dimension of the variable
    int dimension;

    // Boundary flag of the variable
    int bdrflag;

    // List of values of the variable
    float* m_value = nullptr;
  };

  /// \brief Output CFX5 variable object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Variable& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 boundary.
  */
  class gsCFX5Boundary : public gsCFX5Region
  {
  public:
    /// \brief Constructors
    using gsCFX5Region::gsCFX5Region;

    /// \brief Imports CFX5 boundary from CFX5 file
    void importFromCFX5(int id)
    {
      // Check that the boundary ID does not exceed the number of boundaries
      // present in the CFX5 file. Note that if the file has not been
      // initialized then the call to cfxExportBoundaryCount will return 0.
      if (id >= cfxExportBoundaryCount())
        GISMO_ERROR("Boundary ID exceeds number of boundaries in CFX5 file!");

      // Set the boundary ID
      this->id = id;

      // Import the boundary name from CFX5 file
      this->name = cfxExportBoundaryName(id+1);

      // Import the number of nodes from CFX5 file
      this->nodes = cfxExportBoundarySize(id+1, cfxREG_NODES);

      // Import the list of nodes from CFX5 file
      this->m_node_id = cfxCopy(cfxExportBoundaryList(id+1, cfxREG_NODES), this->nodes);

      // Import the number of faces from CFX5 file
      this->faces = cfxExportBoundarySize(id+1, cfxREG_FACES);

      // Import the list of faces from CFX5 file
      this->m_face_id = cfxCopy(cfxExportBoundaryList(id+1, cfxREG_FACES), this->faces);

      // Free internal memory
      cfxExportRegionFree(id+1);
    }

    // Delete inherited functions
    const int getRegion() = delete;
    void setRegion(int id) = delete;

    /// \brief Returns the number of the boundary
    const int getBoundary() const
    { return id; }

    /// \brief Sets the number of the boundary
    void setBoundary(int id)
    { this->id = id; }

    /// \brief Returns hash of the boundary
    size_t getHash() const
    {
      std::string region = util::to_string(getBoundary())+getName()+
        util::to_string(getNodeNumber())+util::to_string(getFaceNumber());

      size_t hash = std::hash<std::string>{}(region);

      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(nodeit->getHash()));

      for (auto faceit = face_cbegin(); faceit != face_cend(); faceit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*faceit).getHash()));

      return hash;
    }

    /// \brief Prints the boundary
    std::ostream& print(std::ostream& os) const
    {
      os << "CFX5Boundary (" << getBoundary() << "): " << getName() << "\n"
         << std::setw(20)    << getNodeNumber()        << "    nodes\n"
         << std::setw(20)    << getFaceNumber()        << "    faces";
      return os;
    }

    /// \brief Prints the boundary as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;
      print(os);

      os << "\n========== NODES IN BOUNDARY ==========\n";
      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        os << (*nodeit).detail() << std::endl;

      os << "\n========== FACES IN BOUNDARY ==========\n";
      for (auto faceit = face_cbegin(); faceit != face_cend(); faceit++)
        os << (*faceit).detail() << std::endl;

      return os.str();
    }

  };

  /// \brief Output CFX5 boundary object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Boundary& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 domain.
  */
  class gsCFX5Domain
  {
  public:
    /// \brief Constructor (default)
    gsCFX5Domain()
      :
      id{0},
      name{"Default Domain"},
      m_node{nullptr},
      m_element{nullptr}
    {
      for (int i=0; i<cfxCNT_SIZE; i++)
        counts[i] = 0;
    }

    /// \brief Constructor (copy)
    gsCFX5Domain(const gsCFX5Domain& other)
      :
      id{other.id},
      name{other.name},
      m_node{cfxCopy(other.getNodes(), other.getNodeNumber())},
      m_element{cfxCopy(other.getElements(), other.getElementNumber())}
    {
#if !defined(NDEBUG)
      gsDebug << "Copying gsCFX5Domain\n";
#endif
      for (int i=0; i<cfxCNT_SIZE; i++)
        counts[i] = other.counts[i];

      std::copy(other.m_region.begin(),
                other.m_region.end(),
                std::back_inserter(m_region));

      std::copy(other.m_volume.begin(),
                other.m_volume.end(),
                std::back_inserter(m_volume));
    }

    /// \brief Constructor (move)
    gsCFX5Domain(gsCFX5Domain&& other)
      :
      id{other.id},
                 name{other.name},
                 m_node{other.m_node},
                 m_element{other.m_element}
    {
#if !defined(NDEBUG)
      gsDebug << "Moving gsCFX5Domain\n";
#endif
      other.id        = 0;
      other.name      = "Default Domain";
      other.m_node    = nullptr;
      other.m_element = nullptr;

      for (int i=0; i<cfxCNT_SIZE; i++)
        {
          counts[i] = other.counts[i];
          other.counts[i] = 0;
        }

      std::swap(m_region, other.m_region); other.m_region.clear();
      std::swap(m_volume, other.m_volume); other.m_volume.clear();
    }

    /// \brief Destructor
    ~gsCFX5Domain()
    {
      id   = 0;
      name = "Default Domain";
      delete[] m_node; m_node = nullptr;
      delete[] m_element; m_element = nullptr;
    }

    /// \brief Copy assignment operator
    gsCFX5Domain& operator=(const gsCFX5Domain& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Copy assigning gsCFX5Domain\n";
#endif
          id   = other.id;
          name = other.name;

          delete[] m_node;    m_node    = cfxCopy(other.getNodes(), other.getNodeNumber());
          delete[] m_element; m_element = cfxCopy(other.getElements(), other.getElementNumber());

          for (int i=0; i<cfxCNT_SIZE; i++)
            counts[i] = other.counts[i];

          m_region.clear();
          std::copy(other.m_region.begin(),
                    other.m_region.end(),
                    std::back_inserter(m_region));

          m_volume.clear();
          std::copy(other.m_volume.begin(),
                    other.m_volume.end(),
                    std::back_inserter(m_volume));
        }
      return *this;
    }

    /// \brief Move assignment operator
    gsCFX5Domain& operator=(gsCFX5Domain&& other)
    {
      if (this != &other)
        {
#if !defined(NDEBUG)
          gsDebug << "Move assigning gsCFX5Domain\n";
#endif
          id   = other.id;   other.id   = 0;
          name = other.name; other.name = "Default Domain";

          delete[] m_node;    m_node    = other.m_node;    other.m_node = nullptr;
          delete[] m_element; m_element = other.m_element; other.m_element = nullptr;

          for (int i=0; i<cfxCNT_SIZE; i++)
            {
              counts[i] = other.counts[i];
              other.counts[i] = 0;
            }

          std::swap(m_region, other.m_region); other.m_region.clear();
          std::swap(m_volume, other.m_volume); other.m_volume.clear();
        }
      return *this;
    }

    /// \brief Returns swapped domain
    gsCFX5Domain& swap(gsCFX5Domain& other)
    {
      std::swap(id,         other.id);
      std::swap(name,       other.name);
      std::swap(counts,     other.counts);
      std::swap(m_node,     other.m_node);
      std::swap(m_element,  other.m_element);
      std::swap(m_region,   other.m_region);
      std::swap(m_boundary, other.m_boundary);
      std::swap(m_volume,   other.m_volume);
      std::swap(m_variable, other.m_variable);
      return *this;
    }

  public:
    /// \brief Allocates memory for the list of nodes and returns a
    /// pointer to it
    gsCFX5Node* createNodes(int nodes)
    {
      delete[] m_node;
      m_node = new gsCFX5Node[nodes];

      return m_node;
    }

    /// \brief Allocates memory for the list of nodes, updates the
    /// internal nodes counter, and returns a pointer to list of nodes
    gsCFX5Node* addNodes(int nodes)
    {
      counts[cfxCNT_NODE] = nodes;
      return createNodes(nodes);
    }

    /// \brief Allocates memory for the list of elements and returns a
    /// pointer to it
    gsCFX5Element* createElements(int elements)
    {
      delete[] m_element;
      m_element = new gsCFX5Element[elements];
      return m_element;
    }

    /// \brief Allocates memory for the list of elements, updates the
    /// internal element counter, and returns a pointer to the list of
    /// element
    gsCFX5Element* addElements(int elements,
                               int tetrahedrals = 0,
                               int pyramids     = 0,
                               int prisms       = 0,
                               int hexahedrals  = 0)
    {
      counts[cfxCNT_ELEMENT] = elements;
      counts[cfxCNT_TET]     = tetrahedrals;
      counts[cfxCNT_PYR]     = pyramids;
      counts[cfxCNT_WDG]     = prisms;
      counts[cfxCNT_HEX]     = hexahedrals;
      return createElements(elements);
    }

  public:
    /// \brief Creates a new region and returns a pointer to it
    gsCFX5Region* createRegion(gsCFX5Domain* m_domain)
    {
      m_region.push_back(std::move(new gsCFX5Region{m_domain}));
      return m_region.back();
    }

    /// \brief Adds a new region, updates the internal region counter,
    /// and returns a pointer to the region object
    gsCFX5Region* addRegion(gsCFX5Domain* m_domain)
    {
      counts[cfxCNT_REGION]++;
      return createRegion(m_domain);
    }

    /// \brief Creates a new boundary region and returns a pointer to
    /// it
    gsCFX5Boundary* createBoundary(gsCFX5Domain* m_domain)
    {
      m_boundary.push_back(std::move(new gsCFX5Boundary{m_domain}));
      return m_boundary.back();
    }

    /// \brief Adds a new boundary region, updates the internal
    /// boundary counter, and returns pointer to the boundary object
    gsCFX5Boundary* addBoundary(gsCFX5Domain* m_domain)
    {
      counts[cfxCNT_BOUNDARY]++;
      return createBoundary(m_domain);
    }

    /// \brief Creates a new volume and returns a pointer to it
    gsCFX5Volume* createVolume(gsCFX5Domain* m_domain)
    {
      m_volume.push_back(std::move(new gsCFX5Volume{m_domain}));
      return m_volume.back();
    }

    /// \brief Adds a new volume, updates the internal volume counter,
    /// and returns a pointer to the volume object
    gsCFX5Volume* addVolume(gsCFX5Domain* m_domain)
    {
      counts[cfxCNT_VOLUME]++;
      return createVolume(m_domain);
    }

    /// \brief Creates a new variable and returns a pointer to it
    gsCFX5Variable* createVariable()
    {
      m_variable.push_back(std::move(new gsCFX5Variable{}));
      return m_variable.back();
    }

    /// \brief Adds a new variable, updates the internal variable
    /// counter, and returns pointer to the variable object
    gsCFX5Variable* addVariable()
    {
      counts[cfxCNT_VARIABLE]++;
      return createVariable();
    }

  public:
    /// \brief Imports CFX5 domain from CFX5 file
    void importFromCFX5(int id)
    {
      // Check that the domain ID does not exceed the number of domains
      // present in the CFX5 file. Note that if the file has not been
      // initialized then the call to cfxExportZoneCount will return 0.
      if (id >= cfxExportZoneCount())
        GISMO_ERROR("Domain ID exceeds number of domains in CFX5 file!");

      // Set the domain ID
      this->id = id;

      // Import the domain name from CFX5 file
      this->name = cfxExportZoneName(id+1);

      // Import the dimensions from CFX5 file
      if (cfxExportZoneSet(id+1, counts) < 0)
        GISMO_ERROR("gsCFX5Domain cannot set the zone " + std::to_string(id+1) + ".");

      // Import nodes from CFX5 file
      m_node = cfxCopy((gsCFX5Node*)cfxExportNodeList(), getNodeNumber());

      // Import elements from CFX5 file
      m_element = cfxCopy((gsCFX5Element*) cfxExportElementList(), getElementNumber());

      // Import regions from CFX5 file
      m_region.reserve(getRegionNumber());
      for (int id=0; id<getRegionNumber(); id++)
        {
          createRegion(this)->importFromCFX5(id);
        }

      // Import boundaries from CFX5 file
      m_boundary.reserve(getBoundaryNumber());
      for (int id=0; id<getBoundaryNumber(); id++)
        {
          createBoundary(this)->importFromCFX5(id);
        }

      // Import volumes from CFX5 file
      m_volume.reserve(getVolumeNumber());
      for (int id=0; id<getVolumeNumber(); id++)
        {
          createVolume(this)->importFromCFX5(id);
        }

      // Import variables from CFX5 file
      m_variable.reserve(getVariableNumber());
      for (int id=0; id<getVariableNumber(); id++)
        {
          createVariable()->importFromCFX5(id);
        }

      // Free internal memory
      cfxExportZoneFree();
    }

    /// \brief Imports CFX5 domain from ICEM ASCII file
    void importFromICEM(int id, std::ifstream& file)
    {
      if (!file.is_open())
        GISMO_ERROR("ICEM file is not ready for import!");

      // Set the domain ID
      this->id = id;

      // Import header from ICEM file
      std::string buffer;
      std::getline(file, buffer); // 1128683573
      std::getline(file, buffer); // Version number: 5.6D

      // Import node, element, region and volume numbers from ICEM file
      int nodes,tetrahedrals,pyramids,prisms,hexahedrals,volumes,regions;
      file >> nodes;
      file >> tetrahedrals;
      file >> prisms;
      file >> hexahedrals;
      file >> pyramids;
      file >> volumes;
      file >> regions;

      // Allocate memory for nodes and elements
      m_node = addNodes(nodes);
      m_element = addElements(tetrahedrals+prisms+hexahedrals+pyramids,
                              tetrahedrals, pyramids, prisms, hexahedrals);

      // Initialize iterators over nodes and elements
      auto nodeit = node_begin();
      auto elemit = element_begin();

      // Import node coordinates from ICEM file
      for (int i=0; i<getNodeNumber(); i++)
        {
          file >> nodeit->getX();
          file >> nodeit->getY();
          file >> nodeit->getZ();
          nodeit++;
        }

      // Import tetrahedral elements from ICEM file
      for (int i=0; i<tetrahedrals; i++)
        {
          elemit->setElementType(gsCFX5ElementType::Tetrahedral);
          file >> elemit->getNodeID(0);
          file >> elemit->getNodeID(1);
          file >> elemit->getNodeID(2);
          elemit++;
        }

      // Import prism elements from ICEM file
      for (int i=0; i<prisms; i++)
        {
          elemit->setElementType(gsCFX5ElementType::Prism);
          file >> elemit->getNodeID(0);
          file >> elemit->getNodeID(1);
          file >> elemit->getNodeID(2);
          file >> elemit->getNodeID(3);
          file >> elemit->getNodeID(4);
          file >> elemit->getNodeID(5);
          elemit++;
        }

      // Import hexahedral elements from ICEM file
      for (int i=0; i<hexahedrals; i++)
        {
          elemit->setElementType(gsCFX5ElementType::Hexahedral);
          file >> elemit->getNodeID(0);
          file >> elemit->getNodeID(1);
          file >> elemit->getNodeID(2);
          file >> elemit->getNodeID(3);
          file >> elemit->getNodeID(4);
          file >> elemit->getNodeID(5);
          file >> elemit->getNodeID(6);
          file >> elemit->getNodeID(7);
          elemit++;
        }

      // Import pyramid elements from ICEM file
      for (int i=0; i<pyramids; i++)
        {
          elemit->setElementType(gsCFX5ElementType::Pyramid);
          file >> elemit->getNodeID(0);
          file >> elemit->getNodeID(1);
          file >> elemit->getNodeID(2);
          file >> elemit->getNodeID(3);
          file >> elemit->getNodeID(4);
          elemit++;
        }

      // Import volumes from ICEM file
      m_region.reserve(volumes);
      for (int id=0; id<volumes; id++)
        addVolume(this)->importFromICEM(id, file);

      // Import regions from ICEM file
      m_region.reserve(regions);
      for (int id=0; id<regions; id++)
        addRegion(this)->importFromICEM(id, file);
    }

    /// \brief Exports CFX5 domain to CFX5 file
    void exportToCFX5() const
    {
      // Check that the CFX5 file is ready
      if (cfxImportStatus() < 0)
        GISMO_ERROR("CFX5 file is not ready for import!");

      // Export nodes to CFX5 file
      ID_t id = 1;
      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++, id++)
        if(cfxImportNode(id, nodeit->getX(), nodeit->getY(), nodeit->getZ())
           != id) GISMO_ERROR("An error occured while importing nodes to CFX5 file!");

      // Export elements to CFX5 file
      id = 1;
      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++, id++)
        if(cfxImportElement(id, (cfxImpElemType_t)elemit->getElementType(),
                            (const ID_t*)elemit->swapExportImport().getNodeIDs())
           != id) GISMO_ERROR("An error occured while importing elements to CFX5 file!");

      // Export volumes to CFX5 file
      for (auto volit = volume_cbegin(); volit != volume_cend(); volit++)
        (*volit)->exportToCFX5();

      // Export regions to CFX5 file
      for (auto regit = region_cbegin(); regit != region_cend(); regit++)
        (*regit)->exportToCFX5();

      // No export of boundaries to CFX5 file
      if (getBoundaryNumber() > 0)
        gsDebug << "Boundaries are not exported to CFX5 file.\n";

      // No export of variables to CFX5 file
      if (getVariableNumber() > 0)
        gsDebug << "Variables are not exported to CFX5 file.\n";
    }

    /// \brief Exports CFX5 domain to ICEM ASCII file
    void exportToICEM(std::ofstream& file) const
    {
      if (!file.is_open())
        GISMO_ERROR("ICEM file is not ready for export!");

      // Export header to ICEM file
      file << "1128683573" << std::endl;
      file << "Version number: 5.6D" << std::endl;

      // Export node, element, region and volume numbers to ICEM file
      file << counts[cfxCNT_NODE]
           << " "
           << counts[cfxCNT_TET]
           << " "
           << counts[cfxCNT_WDG]
           << " "
           << counts[cfxCNT_HEX]
           << " "
           << counts[cfxCNT_PYR]
           << " "
           << counts[cfxCNT_VOLUME]
           << " "
           << counts[cfxCNT_REGION]
           << std::endl;

      // Export node coordinates to ICEM file
      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        file << nodeit->getX()
             << " "
             << nodeit->getY()
             << " "
             << nodeit->getZ()
             << std::endl;

      // Export elements to ICEM file
      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++)
        {
          for (int i=0; i<(int)elemit->getElementType(); i++)
            file << elemit->getNodeID(i) << " ";
          file << std::endl;
        }

      // Export volumes to ICEM file
      for (auto volit = volume_cbegin(); volit != volume_cend(); volit++)
        (*volit)->exportToICEM(file);

      // Export regions to ICEM file
      for (auto regit = region_cbegin(); regit != region_cend(); regit++)
        (*regit)->exportToICEM(file);

      // No export of boundaries to ICEM file
      if (getBoundaryNumber() > 0)
        gsDebug << "Boundaries are not exported to ICEM file.\n";

      // No export of variables to ICEM file
      if (getVariableNumber() > 0)
        gsDebug << "Variables are not exported to ICEM file.\n";
    }

  public:
    /// \brief Returns the number of the domain
    const int getDomain() const
    { return id; }

    /// \brief Sets the number of the domain
    void setDomain(int id)
    { this->id = id; }

    /// \brief Returns the name of the domain
    const std::string& getName() const
    { return name; }

    /// \brief Sets the name of the domain
    void setName(const std::string& name)
    { this->name = name; }

    /// \brief Returns the total number of nodes
    const int getNodeNumber() const
    { return counts[cfxCNT_NODE]; }

    /// \brief Returns the total number of elements
    const int getElementNumber() const
    { return counts[cfxCNT_ELEMENT]; }

    /// \brief Returns the total number of tetrahedra elements
    const int getTetrahedraNumber() const
    { return counts[cfxCNT_TET]; }

    /// \brief Returns the total number of prism elements
    const int getPrismNumber() const
    { return counts[cfxCNT_WDG]; }

    /// \brief Returns the total number of hexahedra elements
    const int getHexahedraNumber() const
    { return counts[cfxCNT_HEX]; }

    /// \brief Returns the total number of pyramid elements
    const int getPyramidNumber() const
    { return counts[cfxCNT_PYR]; }

    /// \brief Returns the  total number of regions
    const int getRegionNumber() const
    { return counts[cfxCNT_REGION]; }

    /// \brief Returns the  total number of boundaries
    const int getBoundaryNumber() const
    { return counts[cfxCNT_BOUNDARY]; }

    /// \brief Returns the total number of volumes
    const int getVolumeNumber() const
    { return counts[cfxCNT_VOLUME]; }

    /// \brief Returns the  total number of variables
    const int getVariableNumber() const
    { return counts[cfxCNT_VARIABLE]; }

    /// \brief Calculates the number of different element types
    void calcElementTypes()
    {
      counts[cfxCNT_TET] = 0;
      counts[cfxCNT_PYR] = 0;
      counts[cfxCNT_WDG] = 0;
      counts[cfxCNT_HEX] = 0;

      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++)
        {
          switch (elemit->getElementType())
            {
            case gsCFX5ElementType::Tetrahedral :
              counts[cfxCNT_TET]++;
              break;
            case gsCFX5ElementType::Pyramid :
              counts[cfxCNT_PYR]++;
              break;
            case gsCFX5ElementType::Prism :
              counts[cfxCNT_WDG]++;
              break;
            case gsCFX5ElementType::Hexahedral :
              counts[cfxCNT_HEX]++;
              break;
            default:
              GISMO_ERROR("Invalid element type!");
            }
        }
    }

  public:
    /// \brief Returns constant pointer to the node array
    const gsCFX5Node* getNodes() const
    {
      return m_node;
    }

    /// \brief Returns pointer to the node array
    gsCFX5Node* getNodes()
    {
      return m_node;
    }

    /// \brief Returns constant reference to the idx-th node
    const gsCFX5Node& getNode(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
      return m_node[idx];
    }

    /// \brief Returns reference to the idx-th node
    gsCFX5Node& getNode(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
      return m_node[idx];
    }

  public:
    /// \brief Returns constant pointer to the element array
    const gsCFX5Element* getElements() const
    {
      return m_element;
    }

    /// \brief Returns pointer to element array
    gsCFX5Element* getElements()
    {
      return m_element;
    }

    /// \brief Returns constant reference to the idx-th element
    const gsCFX5Element& getElement(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getElementNumber()), "Invalid element index.");
      return m_element[idx];
    }

    /// \brief Returns reference to the idx-th element
    gsCFX5Element& getElement(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getElementNumber()), "Invalid element index.");
      return m_element[idx];
    }

    /// \brief Returns hash of the domain
    size_t getHash(bool meshOnly=false) const
    {
      std::string domain = util::to_string(getDomain())+getName()+
        util::to_string(getNodeNumber())+
        util::to_string(getElementNumber())+
        util::to_string(getRegionNumber())+
        util::to_string(getBoundaryNumber())+
        util::to_string(getVolumeNumber());

      if (!meshOnly) domain += util::to_string(getVariableNumber());

      size_t hash = std::hash<std::string>{}(domain);

      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(nodeit->getHash()));

      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string(elemit->getHash()));

      for (auto regit = region_cbegin(); regit != region_cend(); regit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*regit)->getHash()));

      for (auto bdrit = boundary_cbegin(); bdrit != boundary_cend(); bdrit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*bdrit)->getHash()));

      for (auto volit = volume_cbegin(); volit != volume_cend(); volit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*volit)->getHash()));

      if (!meshOnly) {
        for (auto varit = variable_cbegin(); varit != variable_cend(); varit++)
          hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*varit)->getHash()));
      }

      return hash;
    }

    /// \brief Prints the domain
    std::ostream& print(std::ostream& os) const
    {
      os << "CFX5Domain (" << getDomain() << "): "   << getName() << "\n"
         << std::setw(20)  << getNodeNumber()        << "    nodes\n"
         << std::setw(20)  << getElementNumber()     << "    elements\n"
         << std::setw(20)  << getRegionNumber()      << "    regions\n"
         << std::setw(20)  << getBoundaryNumber()    << "    boundaries\n"
         << std::setw(20)  << getVolumeNumber()      << "    volumes\n"
         << std::setw(20)  << getVariableNumber()    << "    variables\n";

      os << "\n========== REGIONS ==========\n";
      for (auto regit = region_cbegin(); regit != region_cend(); regit++)
        os << **regit << std::endl;

      os << "\n========== BOUNDARIES ==========\n";
      for (auto bdrit = boundary_cbegin(); bdrit != boundary_cend(); bdrit++)
        os << **bdrit << std::endl;

      os << "\n========== VOLUMES ==========\n";
      for (auto volit = volume_cbegin(); volit != volume_cend(); volit++)
        os << **volit << std::endl;

      os << "\n========== VARIABLES ==========\n";
      for (auto varit = variable_cbegin(); varit != variable_cend(); varit++)
        os << **varit << std::endl;

      return os;
    }

    /// \brief Prints the domain as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;

      os << "CFX5Domain (" << getDomain() << "): "   << getName() << "\n"
         << std::setw(20)  << getNodeNumber()        << "    nodes\n"
         << std::setw(20)  << getElementNumber()     << "    elements\n"
         << std::setw(20)  << getRegionNumber()      << "    regions\n"
         << std::setw(20)  << getBoundaryNumber()    << "    boundaries\n"
         << std::setw(20)  << getVolumeNumber()      << "    volumes\n"
         << std::setw(20)  << getVariableNumber()    << "    variables\n";

      os << "\n========== NODES IN DOMAIN ==========\n";
      for (auto nodeit = node_cbegin(); nodeit != node_cend(); nodeit++)
        os << (*nodeit).detail() << std::endl;

      os << "\n========== ELEMENTS DOMAIN ==========\n";
      for (auto elemit = element_cbegin(); elemit != element_cend(); elemit++)
        os << (*elemit).detail() << std::endl;

      os << "\n========== REGIONS ==========\n";
      for (auto regit = region_cbegin(); regit != region_cend(); regit++)
        os << (**regit).detail() << std::endl;

      os << "\n========== BOUNDARIES ==========\n";
      for (auto bdrit = boundary_cbegin(); bdrit != boundary_cend(); bdrit++)
        os << (**bdrit).detail() << std::endl;

      os << "\n========== VOLUMES ==========\n";
      for (auto volit = volume_cbegin(); volit != volume_cend(); volit++)
        os << (**volit).detail() << std::endl;

      os << "\n========== VARIABLES ==========\n";
      for (auto varit = variable_cbegin(); varit != variable_cend(); varit++)
        os << (**varit).detail() << std::endl;

      return os.str();
    }

    /// \brief Prints the difference between the domain and the other
    /// domain as a string
    std::string diff(const gsCFX5Domain& other) const
    {
      std::ostringstream os;

      if (getDomain() != other.getDomain())
        os << "CFX5Domain domain\n< "
           << getDomain()
           << "\n---\n"
           << "> "
           << other.getDomain()
           << "\n";

      if (getName() != other.getName())
        os << "CFX5Domain name\n< "
           << getName()
           << "\n---\n"
           << "> "
           << other.getName()
           << "\n";

      if (getNodeNumber() != other.getNodeNumber())
        os << "CFX5Domain nodes\n< "
           << getNodeNumber()
           << "\n---\n"
           << "> "
           << other.getNodeNumber()
           << "\n";

      if (getElementNumber() != other.getElementNumber())
        os << "CFX5Domain elements\n< "
           << getElementNumber()
           << "\n---\n"
           << "> "
           << other.getElementNumber()
           << "\n";

      if (getRegionNumber() != other.getRegionNumber())
        os << "CFX5Domain regions\n< "
           << getRegionNumber()
           << "\n---\n"
           << "> "
           << other.getRegionNumber()
           << "\n";

      if (getBoundaryNumber() != other.getBoundaryNumber())
        os << "CFX5Domain boundarys\n< "
           << getBoundaryNumber()
           << "\n---\n"
           << "> "
           << other.getBoundaryNumber()
           << "\n";

      if (getVolumeNumber() != other.getVolumeNumber())
        os << "CFX5Domain volumes\n< "
           << getVolumeNumber()
           << "\n---\n"
           << "> "
           << other.getVolumeNumber()
           << "\n";

      if (getVariableNumber() != other.getVariableNumber())
        os << "CFX5Domain variables\n< "
           << getVariableNumber()
           << "\n---\n"
           << "> "
           << other.getVariableNumber()
           << "\n";

      os << "\n========== NODES IN DOMAIN ==========\n";
      for (auto nodeit = node_cbegin(), other_nodeit = other.node_cbegin();
           nodeit != node_cend() && other_nodeit != other.node_cend();
           nodeit++, other_nodeit++)
        os << (*nodeit).diff(*other_nodeit);

      os << "\n========== ELEMENTS DOMAIN ==========\n";
      for (auto elemit = element_cbegin(), other_elemit = other.element_cbegin();
           elemit != element_cend() && other_elemit != other.element_cend();
           elemit++, other_elemit++)
        os << (*elemit).diff(*other_elemit);

      os << "\n========== REGIONS ==========\n";
      for (auto regit = region_cbegin(), other_regit = other.region_cbegin();
           regit != region_cend() && other_regit != other.region_cend();
           regit++, other_regit++)
        os << (**regit).diff(**other_regit);


      os << "\n========== BOUNDARIES ==========\n";
      for (auto bdrit = boundary_cbegin(), other_bdrit = other.boundary_cbegin();
           bdrit != boundary_cend() && other_bdrit != other.boundary_cend();
           bdrit++, other_bdrit++)
        os << (**bdrit).diff(**other_bdrit);

      os << "\n========== VOLUMES ==========\n";
      for (auto volit = volume_cbegin(), other_volit = other.volume_cbegin();
           volit != volume_cend() && other_volit != other.volume_cend();
           volit++, other_volit++)
        os << (**volit).diff(**other_volit);

      os << "\n========== VARIABLES ==========\n";
      for (auto varit = variable_cbegin(), other_varit = other.variable_cbegin();
           varit != variable_cend() && other_varit != other.variable_cend();
           varit++, other_varit++)
        os << (**varit).diff(**other_varit);

      return os.str();
    }

  public:
    typedef gsCFX5GenericIterator<gsCFX5Node>             node_iterator;
    typedef gsCFX5GenericConstIterator<gsCFX5Node>        node_const_iterator;
    typedef gsCFX5GenericReverseIterator<gsCFX5Node>      node_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<gsCFX5Node> node_const_reverse_iterator;

    /// Returns an iterator to the beginning of the node list
    node_iterator node_begin()
    {
      return node_iterator(m_node);
    }

    /// Returns a const-iterator to the beginning of the node list
    node_const_iterator node_cbegin() const
    {
      return node_const_iterator(m_node);
    }

    /// Returns a reverse-iterator to the beginning of the node list
    node_reverse_iterator node_rbegin()
    {
      return node_reverse_iterator(m_node + getNodeNumber());
    }

    /// Returns a const-reverse-iterator to the beginning of the node list
    node_const_reverse_iterator node_crbegin() const
    {
      return node_const_reverse_iterator(m_node + getNodeNumber());
    }

    /// Returns an iterator to the end of the node list
    node_iterator node_end()
    {
      return node_iterator(m_node + getNodeNumber());
    }

    /// Returns a const-iterator to the end of the node list
    node_const_iterator node_cend() const
    {
      return node_const_iterator(m_node + getNodeNumber());
    }

    /// Returns a reverse-iterator to the end of the node list
    node_reverse_iterator node_rend()
    {
      return node_reverse_iterator(m_node);
    }

    /// Returns a const-reverse-iterator to the end of the node list
    node_const_reverse_iterator node_crend() const
    {
      return node_const_reverse_iterator(m_node);
    }

  public:
    typedef gsCFX5GenericIterator<gsCFX5Element>             element_iterator;
    typedef gsCFX5GenericConstIterator<gsCFX5Element>        element_const_iterator;
    typedef gsCFX5GenericReverseIterator<gsCFX5Element>      element_reverse_iterator;
    typedef gsCFX5GenericConstReverseIterator<gsCFX5Element> element_const_reverse_iterator;

    /// Returns an iterator to the beginning of the element list
    element_iterator element_begin()
    {
      return element_iterator(m_element);
    }

    /// Returns a const-iterator to the beginning of the element list
    element_const_iterator element_cbegin() const
    {
      return element_const_iterator(m_element);
    }

    /// Returns a reverse-iterator to the beginning of the element list
    element_reverse_iterator element_rbegin()
    {
      return element_reverse_iterator(m_element + getElementNumber());
    }

    /// Returns a const-reverse-iterator to the beginning of the element list
    element_const_reverse_iterator element_crbegin() const
    {
      return element_const_reverse_iterator(m_element + getElementNumber());
    }

    /// Returns an iterator to the end of the element list
    element_iterator element_end()
    {
      return element_iterator(m_element + getElementNumber());
    }

    /// Returns a const-iterator to the end of the element list
    element_const_iterator element_cend() const
    {
      return element_const_iterator(m_element + getElementNumber());
    }

    /// Returns a reverse-iterator to the end of the element list
    element_reverse_iterator element_rend()
    {
      return element_reverse_iterator(m_element);
    }

    /// Returns a const-reverse-iterator to the end of the element list
    element_const_reverse_iterator element_crend() const
    {
      return element_const_reverse_iterator(m_element);
    }

  public:
    typedef typename std::vector<gsCFX5Region*>     Region;
    typedef typename Region::iterator               region_iterator;
    typedef typename Region::const_iterator         region_const_iterator;
    typedef typename Region::reverse_iterator       region_reverse_iterator;
    typedef typename Region::const_reverse_iterator region_const_reverse_iterator;

    /// Returns an iterator to the beginning of the region vector
    region_iterator region_begin()
    { return m_region.begin(); }

    /// Returns a const-iterator to the beginning of the region vector
    region_const_iterator region_cbegin() const
    { return m_region.cbegin(); }

    /// Returns a reverse-iterator to the beginning of the region vector
    region_reverse_iterator region_rbegin()
    { return m_region.rbegin(); }

    /// Returns a const-reverse-iterator to the beginning of the region vector
    region_const_reverse_iterator region_crbegin() const
    { return m_region.crbegin(); }

    /// Returns an iterator to the end of the region vector
    region_iterator region_end()
    { return m_region.end(); }

    /// Returns a const-iterator to the end of the region vector
    region_const_iterator region_cend() const
    { return m_region.cend(); }

    /// Returns a reverse-iterator to the end of the region vector
    region_reverse_iterator region_rend()
    { return m_region.rend(); }

    /// Returns a const-reverse-iterator to the end of the region vector
    region_const_reverse_iterator region_crend() const
    { return m_region.crend(); }

  public:
    /// Returns constant reference to region vector
    const Region& getRegion() const
    { return m_region; }

    /// Returns reference to region vector
    Region& getRegion()
    { return m_region; }

    /// Returns constant reference to the idx-th region
    const gsCFX5Region* getRegion(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getRegionNumber()), "Invalid region index.");
      return m_region[idx];
    }

    /// Returns reference to the idx-th region
    gsCFX5Region* getRegion(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getRegionNumber()), "Invalid region index.");
      return m_region[idx];
    }

    /// Returns constant pointer to region with given ID
    /// If the region ID is not found a nullptr is returned
    const gsCFX5Region* searchRegion(int region) const
    {
      for (auto it = region_cbegin(); it != region_cend(); it++)
        if ((*it)->getRegion() == region)
          return *it;
      return nullptr;
    }

    /// Returns pointer to region with given ID
    /// If the region ID is not found a nullptr is returned
    gsCFX5Region* searchRegion(int region)
    {
      for (auto it = region_begin(); it != region_end(); it++)
        if ((*it)->getRegion() == region)
          return *it;
      return nullptr;
    }

    /// Returns constant pointer to region with given name
    /// If the region name is not found a nullptr is returned
    const gsCFX5Region* searchRegion(const std::string& name) const
    {
      for (auto it = region_cbegin(); it != region_cend(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

    /// Returns pointer to region with given name
    /// If the region name is not found a nullptr is returned
    gsCFX5Region* searchRegion(const std::string& name)
    {
      for (auto it = region_begin(); it != region_end(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

  public:
    typedef typename std::vector<gsCFX5Boundary*>     Boundary;
    typedef typename Boundary::iterator               boundary_iterator;
    typedef typename Boundary::const_iterator         boundary_const_iterator;
    typedef typename Boundary::reverse_iterator       boundary_reverse_iterator;
    typedef typename Boundary::const_reverse_iterator boundary_const_reverse_iterator;

    /// Returns an iterator to the beginning of the boundary vector
    boundary_iterator boundary_begin()
    { return m_boundary.begin(); }

    /// Returns a const-iterator to the beginning of the boundary vector
    boundary_const_iterator boundary_cbegin() const
    { return m_boundary.cbegin(); }

    /// Returns a reverse-iterator to the beginning of the boundary vector
    boundary_reverse_iterator boundary_rbegin()
    { return m_boundary.rbegin(); }

    /// Returns a const-reverse-iterator to the beginning of the boundary vector
    boundary_const_reverse_iterator boundary_crbegin() const
    { return m_boundary.crbegin(); }

    /// Returns an iterator to the end of the boundary vector
    boundary_iterator boundary_end()
    { return m_boundary.end(); }

    /// Returns a const-iterator to the end of the boundary vector
    boundary_const_iterator boundary_cend() const
    { return m_boundary.cend(); }

    /// Returns a reverse-iterator to the end of the boundary vector
    boundary_reverse_iterator boundary_rend()
    { return m_boundary.rend(); }

    /// Returns a const-reverse-iterator to the end of the boundary vector
    boundary_const_reverse_iterator boundary_crend() const
    { return m_boundary.crend(); }

  public:
    /// Returns constant reference to boundary vector
    const Boundary& getBoundary() const
    { return m_boundary; }

    /// Returns reference to boundary vector
    Boundary& getBoundary()
    { return m_boundary; }

    /// Returns constant reference to the idx-th boundary
    const gsCFX5Boundary* getBoundary(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getBoundaryNumber()), "Invalid boundary index.");
      return m_boundary[idx];
    }

    /// Returns reference to the idx-th boundary
    gsCFX5Boundary* getBoundary(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getBoundaryNumber()), "Invalid boundary index.");
      return m_boundary[idx];
    }

    /// Returns constant pointer to boundary with given ID
    /// If the boundary ID is not found a nullptr is returned
    const gsCFX5Boundary* searchBoundary(int boundary) const
    {
      for (auto it = boundary_cbegin(); it != boundary_cend(); it++)
        if ((*it)->getBoundary() == boundary)
          return *it;
      return nullptr;
    }

    /// Returns pointer to boundary with given ID
    /// If the boundary ID is not found a nullptr is returned
    gsCFX5Boundary* searchBoundary(int boundary)
    {
      for (auto it = boundary_begin(); it != boundary_end(); it++)
        if ((*it)->getBoundary() == boundary)
          return *it;
      return nullptr;
    }

    /// Returns constant pointer to boundary with given name
    /// If the boundary name is not found a nullptr is returned
    const gsCFX5Boundary* searchBoundary(const std::string& name) const
    {
      for (auto it = boundary_cbegin(); it != boundary_cend(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

    /// Returns pointer to boundary with given name
    /// If the boundary name is not found a nullptr is returned
    gsCFX5Boundary* searchBoundary(const std::string& name)
    {
      for (auto it = boundary_begin(); it != boundary_end(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

  public:
    typedef typename std::vector<gsCFX5Volume*>     Volume;
    typedef typename Volume::iterator               volume_iterator;
    typedef typename Volume::const_iterator         volume_const_iterator;
    typedef typename Volume::reverse_iterator       volume_reverse_iterator;
    typedef typename Volume::const_reverse_iterator volume_const_reverse_iterator;

    /// Returns an iterator to the beginning of the volume vector
    volume_iterator volume_begin()
    { return m_volume.begin(); }

    /// Returns a const-iterator to the beginning of the volume vector
    volume_const_iterator volume_cbegin() const
    { return m_volume.cbegin(); }

    /// Returns a reverse-iterator to the beginning of the volume vector
    volume_reverse_iterator volume_rbegin()
    { return m_volume.rbegin(); }

    /// Returns a const-reverse-iterator to the beginning of the volume vector
    volume_const_reverse_iterator volume_crbegin() const
    { return m_volume.crbegin(); }

    /// Returns an iterator to the end of the volume vector
    volume_iterator volume_end()
    { return m_volume.end(); }

    /// Returns a const-iterator to the end of the volume vector
    volume_const_iterator volume_cend() const
    { return m_volume.cend(); }

    /// Returns a reverse-iterator to the end of the volume vector
    volume_reverse_iterator volume_rend()
    { return m_volume.rend(); }

    /// Returns a const-reverse-iterator to the end of the volume vector
    volume_const_reverse_iterator volume_crend() const
    { return m_volume.crend(); }

  public:
    /// Returns constant reference to volume vector
    const Volume& getVolume() const
    { return m_volume; }

    /// Returns reference to volume vector
    Volume& getVolume()
    { return m_volume; }

    /// Returns constant reference to the idx-th volume
    const gsCFX5Volume* getVolume(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getVolumeNumber()), "Invalid volume index.");
      return m_volume[idx];
    }

    /// Returns reference to the idx-th volume
    gsCFX5Volume* getVolume(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getVolumeNumber()), "Invalid volume index.");
      return m_volume[idx];
    }

    /// Returns constant pointer to volume with given ID
    /// If the volume ID is not found a nullptr is returned
    const gsCFX5Volume* searchVolume(int id) const
    {
      for (auto it = volume_cbegin(); it != volume_cend(); it++)
        if ((*it)->getVolume() == id)
          return *it;
      return nullptr;
    }

    /// Returns pointer to volume with given ID
    /// If the volume ID is not found a nullptr is returned
    gsCFX5Volume* searchVolume(int id)
    {
      for (auto it = volume_begin(); it != volume_end(); it++)
        if ((*it)->getVolume() == id)
          return *it;
      return nullptr;
    }

    /// Returns constant pointer to volume with given name
    /// If the volume name is not found a nullptr is returned
    const gsCFX5Volume* searchVolume(const std::string& name) const
    {
      for (auto it = volume_cbegin(); it != volume_cend(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

    /// Returns pointer to volume with given name
    /// If the volume name is not found a nullptr is returned
    gsCFX5Volume* searchVolume(const std::string& name)
    {
      for (auto it = volume_begin(); it != volume_end(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

  public:
    typedef typename std::vector<gsCFX5Variable*>     Variable;
    typedef typename Variable::iterator               variable_iterator;
    typedef typename Variable::const_iterator         variable_const_iterator;
    typedef typename Variable::reverse_iterator       variable_reverse_iterator;
    typedef typename Variable::const_reverse_iterator variable_const_reverse_iterator;

    /// Returns an iterator to the beginning of the variable vector
    variable_iterator variable_begin()
    { return m_variable.begin(); }

    /// Returns a const-iterator to the beginning of the variable vector
    variable_const_iterator variable_cbegin() const
    { return m_variable.cbegin(); }

    /// Returns a reverse-iterator to the beginning of the variable vector
    variable_reverse_iterator variable_rbegin()
    { return m_variable.rbegin(); }

    /// Returns a const-reverse-iterator to the beginning of the variable vector
    variable_const_reverse_iterator variable_crbegin() const
    { return m_variable.crbegin(); }

    /// Returns an iterator to the end of the variable vector
    variable_iterator variable_end()
    { return m_variable.end(); }

    /// Returns a const-iterator to the end of the variable vector
    variable_const_iterator variable_cend() const
    { return m_variable.cend(); }

    /// Returns a reverse-iterator to the end of the variable vector
    variable_reverse_iterator variable_rend()
    { return m_variable.rend(); }

    /// Returns a const-reverse-iterator to the end of the variable vector
    variable_const_reverse_iterator variable_crend() const
    { return m_variable.crend(); }

  public:
    /// Returns constant reference to variable vector
    const Variable& getVariable() const
    { return m_variable; }

    /// Returns reference to variable vector
    Variable& getVariable()
    { return m_variable; }

    /// Returns constant reference to the idx-th variable
    const gsCFX5Variable* getVariable(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getVariableNumber()), "Invalid variable index.");
      return m_variable[idx];
    }

    /// Returns reference to the idx-th variable
    gsCFX5Variable* getVariable(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getVariableNumber()), "Invalid variable index.");
      return m_variable[idx];
    }

    /// Returns constant pointer to variable with given ID
    /// If the variable ID is not found a nullptr is returned
    const gsCFX5Variable* searchVariable(int variable) const
    {
      for (auto it = variable_cbegin(); it != variable_cend(); it++)
        if ((*it)->getVariable() == variable)
          return *it;
      return nullptr;
    }

    /// Returns pointer to variable with given ID
    /// If the variable ID is not found a nullptr is returned
    gsCFX5Variable* searchVariable(int variable)
    {
      for (auto it = variable_begin(); it != variable_end(); it++)
        if ((*it)->getVariable() == variable)
          return *it;
      return nullptr;
    }

    /// Returns constant pointer to variable with given name
    /// If the variable name is not found a nullptr is returned
    const gsCFX5Variable* searchVariable(const std::string& name, bool alias=false) const
    {
      if (alias)
        {
          for (auto it = variable_cbegin(); it != variable_cend(); it++)
            if ((*it)->getAlias().compare(name) == 0)
              return *it;
        }
      else
        {
          for (auto it = variable_cbegin(); it != variable_cend(); it++)
            if ((*it)->getName().compare(name) == 0)
              return *it;
        }
      return nullptr;
    }

    /// Returns pointer to variable with given name
    /// If the variable name is not found a nullptr is returned
    gsCFX5Variable* searchVariable(const std::string& name, bool alias=false)
    {
      if (alias)
        {
          for (auto it = variable_begin(); it != variable_end(); it++)
            if ((*it)->getAlias().compare(name) == 0)
              return *it;
        }
      else
        {
          for (auto it = variable_begin(); it != variable_end(); it++)
            if ((*it)->getName().compare(name) == 0)
              return *it;
        }
      return nullptr;
    }

  private:
    /// Domain ID
    int id;

    /// Domain name
    std::string name;

    /// List of nodes in the domain
    gsCFX5Node* m_node = nullptr;

    /// List of elements in the domain
    gsCFX5Element* m_element = nullptr;

    /// Element/node counters for the domain
    int counts[cfxCNT_SIZE];

    /// Region vector
    Region m_region;

    /// Boundary vector
    Boundary m_boundary;

    /// Volume vector
    Volume m_volume;

    /// Variable vector
    Variable m_variable;
  };

  /// \brief Output CFX5 domain object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5Domain& obj)
  {
    return obj.print(os);
  }

  /**
     \brief A CFX5 result.
  */
  class gsCFX5
  {
  public:

    /// \brief Constructor (default)
    gsCFX5() = default;

    /// \brief Destructor
    ~gsCFX5()
    {
      m_domain.clear();
    };

    /// \brief Returns hash of the CFX5 result
    size_t getHash() const
    {
      size_t hash = std::hash<std::string>{}(util::to_string(getDomainNumber()));

      for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
        hash = std::hash<std::string>{}(util::to_string(hash)+util::to_string((*domit)->getHash()));

      return hash;
    }

  public:
    /// \brief Creates a new domain and returns pointer to it
    gsCFX5Domain* createDomain()
    {
      m_domain.push_back(std::move(new gsCFX5Domain{}));
      return m_domain.back();
    }

  public:
    /// \brief Imports CFX5 results from CFX5 file
    void importFromCFX5(const std::string& filename)
    {
      gsDebug << "Importing CFX5 results from CFX5 file <" << filename << ">\n";
      int domains = cfxExportInit(const_cast<char*>(filename.c_str()), NULL);

      m_domain.reserve(domains);
      for (int domain=0; domain<domains; domain++)
        {
          createDomain()->importFromCFX5(domain);
        }

      cfxExportDone();
    }

    /// \brief Imports CFX5 mesh from ICEM ASCII file
    void importFromICEM(const std::string& filename)
    {
      gsDebug << "Importing CFX5 mesh from ICEM ASCII file <" << filename << ">\n";

      std::ifstream file;
      file.open(filename);
      createDomain()->importFromICEM(0, file);
      file.close();
    }

    /// \brief Exports CFX5 results to CFX5 instance
    void exportToCFX5() const
    {
      gsDebug << "Exporting CFX5 results to CFX5 instance\n";
      cfxImportInit();

      for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
        (*domit)->exportToCFX5();

      cfxImportError(NULL);
      cfxImportDone();

#if !defined(NDEBUG)
      size_t counts[cfxImpCNT_SIZE];
      long bytes = cfxImportTotals(counts);
      gsDebug << "Details of CFX5 export:\n"
              << std::setw(20) << counts[cfxImpCNT_REGION]  << " regions\n"
              << std::setw(20) << counts[cfxImpCNT_NODE]    << " nodes\n"
              << std::setw(20) << counts[cfxImpCNT_UNUSED]  << " unused nodes\n"
              << std::setw(20) << counts[cfxImpCNT_DUP]     << " duplicate nodes\n"
              << std::setw(20) << counts[cfxImpCNT_ELEMENT] << " elements\n"
              << std::setw(20) << counts[cfxImpCNT_TET]     << " tetrahedral elements\n"
              << std::setw(20) << counts[cfxImpCNT_PYR]     << " pyramid elements\n"
              << std::setw(20) << counts[cfxImpCNT_WDG]     << " prism elements\n"
              << std::setw(20) << counts[cfxImpCNT_HEX]     << " hexahedral\n"
              << std::setw(20) << bytes                     << " total bytes\n";
#endif
    }

    /// \brief Exports CFX5 results to CFX5 file
    void exportToCFX5(const std::string& filename) const
    {
      gsDebug << "Exporting CFX5 results to CFX5 file <" << filename << ">\n";
      cfxImportTest(const_cast<char*>(filename.c_str()));

      for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
        (*domit)->exportToCFX5();

      cfxImportError(NULL);
      cfxImportDone();

#if !defined(NDEBUG)
      size_t counts[cfxImpCNT_SIZE];
      long bytes = cfxImportTotals(counts);
      gsDebug << "Details of CFX5 export:\n"
              << std::setw(20) << counts[cfxImpCNT_REGION]  << " regions\n"
              << std::setw(20) << counts[cfxImpCNT_NODE]    << " nodes\n"
              << std::setw(20) << counts[cfxImpCNT_UNUSED]  << " unused nodes\n"
              << std::setw(20) << counts[cfxImpCNT_DUP]     << " duplicate nodes\n"
              << std::setw(20) << counts[cfxImpCNT_ELEMENT] << " elements\n"
              << std::setw(20) << counts[cfxImpCNT_TET]     << " tetrahedral elements\n"
              << std::setw(20) << counts[cfxImpCNT_PYR]     << " pyramid elements\n"
              << std::setw(20) << counts[cfxImpCNT_WDG]     << " prism elements\n"
              << std::setw(20) << counts[cfxImpCNT_HEX]     << " hexahedral\n"
              << std::setw(20) << bytes                     << " total bytes\n";
#endif
    }

    /// \brief Exports CFX5 mesh to ICEM ASCII file
    void exportToICEM(const std::string& filename) const
    {
      gsDebug << "Exporting CFX5 results to CFX5 file <" << filename << ">\n";

      std::ofstream file;
      file.open(filename);
      for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
        (*domit)->exportToICEM(file);
      file.close();
    }

  public:
    typedef typename std::vector<gsCFX5Domain*>     Domain;
    typedef typename Domain::iterator               domain_iterator;
    typedef typename Domain::const_iterator         domain_const_iterator;
    typedef typename Domain::reverse_iterator       domain_reverse_iterator;
    typedef typename Domain::const_reverse_iterator domain_const_reverse_iterator;

    /// Returns an iterator to the beginning of the domain vector
    domain_iterator domain_begin()
    { return m_domain.begin(); }

    /// Returns a const-iterator to the beginning of the domain vector
    domain_const_iterator domain_cbegin() const
    { return m_domain.cbegin(); }

    /// Returns a reverse-iterator to the beginning of the domain vector
    domain_reverse_iterator domain_rbegin()
    { return m_domain.rbegin(); }

    /// Returns a const-reverse-iterator to the beginning of the domain vector
    domain_const_reverse_iterator domain_crbegin() const
    { return m_domain.crbegin(); }

    /// Returns an iterator to the end of the domain vector
    domain_iterator domain_end()
    { return m_domain.end(); }

    /// Returns a const-iterator to the end of the domain vector
    domain_const_iterator domain_cend() const
    { return m_domain.cend(); }

    /// Returns a reverse-iterator to the end of the domain vector
    domain_reverse_iterator domain_rend()
    { return m_domain.rend(); }

    /// Returns a const-reverse-iterator to the end of the domain vector
    domain_const_reverse_iterator domain_crend() const
    { return m_domain.crend(); }

  public:
    /// Returns number of domains
    const int getDomainNumber() const
    { return m_domain.size(); }

    /// Returns constant reference to domain vector
    const Domain& getDomain() const
    { return m_domain; }

    /// Returns reference to domain vector
    Domain& getDomain()
    { return m_domain; }

    /// Returns constant reference to the idx-th domain
    const gsCFX5Domain* getDomain(int idx) const
    { return m_domain[idx]; }

    /// Returns reference to the idx-th domain
    gsCFX5Domain* getDomain(int idx)
    { return m_domain[idx]; }

    /// Returns constant pointer to domain with given ID
    /// If the domain ID is not found a nullptr is returned
    const gsCFX5Domain* searchDomain(int domain) const
    {
      for (auto it = domain_cbegin(); it != domain_cend(); it++)
        if ((*it)->getDomain() == domain)
          return *it;
      return nullptr;
    }

    /// Returns pointer to domain with given ID
    /// If the domain ID is not found a nullptr is returned
    gsCFX5Domain* searchDomain(int domain)
    {
      for (auto it = domain_begin(); it != domain_end(); it++)
        if ((*it)->getDomain() == domain)
          return *it;
      return nullptr;
    }

    /// Returns constant pointer to domain with given name
    /// If the domain name is not found a nullptr is returned
    const gsCFX5Domain* searchDomain(const std::string& name) const
    {
      for (auto it = domain_cbegin(); it != domain_cend(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

    /// Returns pointer to domain with given name
    /// If the domain name is not found a nullptr is returned
    gsCFX5Domain* searchDomain(const std::string& name)
    {
      for (auto it = domain_begin(); it != domain_end(); it++)
        if ((*it)->getName().compare(name) == 0)
          return *it;
      return nullptr;
    }

    /// \brief Prints the CFX5 result object
    std::ostream& print( std::ostream& os ) const
    {
      os << "CFX5 (" << getDomainNumber() << "):\n";

      for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
        os << **domit;

      return os;
    }

    /// \brief Prints the CFX5 result object as a string with extended details
    std::string detail() const
    {
      std::ostringstream os;

      os << "CFX5 (" << getDomainNumber() << "):\n";

      for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
        os << (**domit).detail();

      return os.str();
    }

    /// \brief Prints the difference between the CFX5 result object
    /// and the other CFX result object  as a string
    std::string diff(const gsCFX5& other) const
    {
      std::ostringstream os;

      if (getDomainNumber() != other.getDomainNumber())
        os << "CFX5 domains\n< "
           << getDomainNumber()
           << "\n---\n"
           << "> "
           << other.getDomainNumber()
           << "\n";

      for (auto domit = domain_cbegin(), other_domit = other.domain_cbegin() ;
           domit != domain_cend(), other_domit != other.domain_cend() ;
           domit++, other_domit++)
        os << (**domit).diff(**other_domit);

      return os.str();
    }

  public:

    /// Returns swapped CFX5 result object
    gsCFX5& swap(gsCFX5& other)
    {
      std::swap(m_domain, other.m_domain);
      return *this;
    }

  private:
    /// Domain vector
    Domain m_domain;
  };

  /// \brief Output CFX5 results object to output stream
  std::ostream& operator<<(std::ostream& os, const gsCFX5& obj)
  {
    return obj.print(os);
  }

  /// \brief Returns constant reference to the idx-th node in the region
  const gsCFX5Node& gsCFX5Region::getNode(int idx) const
  {
    GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
    return m_domain->getNode(m_node_id[idx]-1);
  }

  /// \brief Returns reference to the idx-th node in the region
  gsCFX5Node& gsCFX5Region::getNode(int idx)
  {
    GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
    return m_domain->getNode(m_node_id[idx]-1);
  }

  /// \brief Returns constant reference to the element associated with
  /// the idx-th face in the region
  const gsCFX5Element& gsCFX5Region::getElement(int idx) const
  {
    GISMO_ASSERT( (idx >= 0) && (idx < getFaceNumber()), "Invalid face index.");
    return m_domain->getElement(cfxELEMNUM(m_face_id[idx])-1);
  }

  /// \brief Returns reference to the element associated with the
  /// idx-th face in the region
  gsCFX5Element& gsCFX5Region::getElement(int idx)
  {
    GISMO_ASSERT( (idx >= 0) && (idx < getFaceNumber()), "Invalid face index.");
    return m_domain->getElement(cfxELEMNUM(m_face_id[idx])-1);
  }

      /// Returns an iterator to the beginning of the node list
  gsCFX5Region::node_iterator gsCFX5Region::node_begin()
    {
      return gsCFX5Region::node_iterator(m_node_id, m_domain->getNodes());
    }

    /// Returns a const-iterator to the beginning of the node list
  gsCFX5Region::node_const_iterator gsCFX5Region::node_cbegin() const
    {
      return gsCFX5Region::node_const_iterator(m_node_id, m_domain->getNodes());
    }

    /// Returns a reverse-iterator to the beginning of the node list
    gsCFX5Region::node_reverse_iterator gsCFX5Region::node_rbegin()
    {
      return gsCFX5Region::node_reverse_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
    }

    /// Returns a const-reverse-iterator to the beginning of the node list
    gsCFX5Region::node_const_reverse_iterator gsCFX5Region::node_crbegin() const
    {
      return gsCFX5Region::node_const_reverse_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
    }

    /// Returns an iterator to the end of the node list
    gsCFX5Region::node_iterator gsCFX5Region::node_end()
    {
      return gsCFX5Region::node_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
    }

    /// Returns a const-iterator to the end of the node list
    gsCFX5Region::node_const_iterator gsCFX5Region::node_cend() const
    {
      return gsCFX5Region::node_const_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
    }

    /// Returns a reverse-iterator to the end of the node list
    gsCFX5Region::node_reverse_iterator gsCFX5Region::node_rend()
    {
      return gsCFX5Region::node_reverse_iterator(m_node_id, m_domain->getNodes());
    }

    /// Returns a const-reverse-iterator to the end of the node list
    gsCFX5Region::node_const_reverse_iterator gsCFX5Region::node_crend() const
    {
      return gsCFX5Region::node_const_reverse_iterator(m_node_id, m_domain->getNodes());
    }

  /// \brief Dereference operator
  const gsCFX5FaceInRegionConstIterator::value_type
  gsCFX5FaceInRegionConstIterator::operator*() const
  {
    return ptr_region->getFace(ptr);
  }

  /// \brief Dereference operator
  const gsCFX5FaceInRegionConstReverseIterator::value_type
  gsCFX5FaceInRegionConstReverseIterator::operator*() const
  { return ptr_region->getFace(ptr); }

    /// \brief Returns constant reference to the idx-th node
  const gsCFX5Node& gsCFX5Volume::getNode(int idx) const
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
      return m_domain->getNode(m_node_id[idx]-1);
    }

    /// \brief Returns reference to the idx-th node
  gsCFX5Node& gsCFX5Volume::getNode(int idx)
    {
      GISMO_ASSERT( (idx >= 0) && (idx < getNodeNumber()), "Invalid node index.");
      return m_domain->getNode(m_node_id[idx]-1);
    }

  /// Returns an iterator to the beginning of the node list
  gsCFX5Volume::node_iterator gsCFX5Volume::node_begin()
  {
    return gsCFX5Volume::node_iterator(m_node_id, m_domain->getNodes());
  }

  /// Returns a const-iterator to the beginning of the node list
  gsCFX5Volume::node_const_iterator gsCFX5Volume::node_cbegin() const
  {
    return gsCFX5Volume::node_const_iterator(m_node_id, m_domain->getNodes());
  }

  /// Returns a reverse-iterator to the beginning of the node list
  gsCFX5Volume::node_reverse_iterator gsCFX5Volume::node_rbegin()
  {
    return gsCFX5Volume::node_reverse_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
  }

  /// Returns a const-reverse-iterator to the beginning of the node list
  gsCFX5Volume::node_const_reverse_iterator gsCFX5Volume::node_crbegin() const
  {
    return gsCFX5Volume::node_const_reverse_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
  }

  /// Returns an iterator to the end of the node list
  gsCFX5Volume::node_iterator gsCFX5Volume::node_end()
  {
    return gsCFX5Volume::node_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
  }

  /// Returns a const-iterator to the end of the node list
  gsCFX5Volume::node_const_iterator gsCFX5Volume::node_cend() const
  {
    return gsCFX5Volume::node_const_iterator(m_node_id + getNodeNumber(), m_domain->getNodes());
  }

  /// Returns a reverse-iterator to the end of the node list
  gsCFX5Volume::node_reverse_iterator gsCFX5Volume::node_rend()
  {
    return gsCFX5Volume::node_reverse_iterator(m_node_id, m_domain->getNodes());
  }

  /// Returns a const-reverse-iterator to the end of the node list
  gsCFX5Volume::node_const_reverse_iterator gsCFX5Volume::node_crend() const
  {
    return gsCFX5Volume::node_const_reverse_iterator(m_node_id, m_domain->getNodes());
  }

  /// Returns an iterator to the beginning of the element list
  gsCFX5Volume::element_iterator gsCFX5Volume::element_begin()
  {
    return gsCFX5Volume::element_iterator(m_element_id, m_domain->getElements());
  }

  /// Returns a const-iterator to the beginning of the element list
  gsCFX5Volume::element_const_iterator gsCFX5Volume::element_cbegin() const
  {
    return gsCFX5Volume::element_const_iterator(m_element_id, m_domain->getElements());
  }

  /// Returns a reverse-iterator to the beginning of the element list
  gsCFX5Volume::element_reverse_iterator gsCFX5Volume::element_rbegin()
  {
    return gsCFX5Volume::element_reverse_iterator(m_element_id + getElementNumber(), m_domain->getElements());
  }

  /// Returns a const-reverse-iterator to the beginning of the element list
  gsCFX5Volume::element_const_reverse_iterator gsCFX5Volume::element_crbegin() const
  {
    return gsCFX5Volume::element_const_reverse_iterator(m_element_id + getElementNumber(), m_domain->getElements());
  }

  /// Returns an iterator to the end of the element list
  gsCFX5Volume::element_iterator gsCFX5Volume::element_end()
  {
    return gsCFX5Volume::element_iterator(m_element_id + getElementNumber(), m_domain->getElements());
  }

  /// Returns a const-iterator to the end of the element list
  gsCFX5Volume::element_const_iterator gsCFX5Volume::element_cend() const
  {
    return gsCFX5Volume::element_const_iterator(m_element_id + getElementNumber(), m_domain->getElements());
  }

  /// Returns a reverse-iterator to the end of the element list
  gsCFX5Volume::element_reverse_iterator gsCFX5Volume::element_rend()
  {
    return gsCFX5Volume::element_reverse_iterator(m_element_id, m_domain->getElements());
  }

  /// Returns a const-reverse-iterator to the end of the element list
  gsCFX5Volume::element_const_reverse_iterator gsCFX5Volume::element_crend() const
  {
    return gsCFX5Volume::element_const_reverse_iterator(m_element_id, m_domain->getElements());
  }

  /// Convert face ID from export to import format
  void gsCFX5Region::convertFaceIDs()
  {
    for (auto faceidit = faceid_begin(); faceidit != faceid_end(); faceidit++)
      {
        int elem = cfxELEMNUM( *faceidit );
        int face = cfxFACENUM( *faceidit );

        switch (m_domain->getElement(elem-1).getElementType())
          {
          case gsCFX5ElementType::Tetrahedral :
            // Do nothing
            break;

          case gsCFX5ElementType::Pyramid :
            switch (face)
              {
              case 1 :
                *faceidit = cfxFACEID(elem, 5);
                break;
              case 2 :
                *faceidit = cfxFACEID(elem, 3);
                break;
              case 3 :
                *faceidit = cfxFACEID(elem, 2);
                break;
              case 4 :
                // Do nothing
                break;
              case 5 :
                *faceidit = cfxFACEID(elem, 1);
                break;
              default :
                GISMO_ERROR("Invalid face number!");
              }
            break;

          case gsCFX5ElementType::Prism :
            switch (face)
              {
              case 1 :
                *faceidit = cfxFACEID(elem, 4);
                break;
              case 2 :
                *faceidit = cfxFACEID(elem, 3);
                break;
              case 3 :
                *faceidit = cfxFACEID(elem, 5);
                break;
              case 4 :
                *faceidit = cfxFACEID(elem, 1);
                break;
              case 5 :
                *faceidit = cfxFACEID(elem, 2);
                break;
              default :
                GISMO_ERROR("Invalid face number!");
              }
            break;

          case gsCFX5ElementType::Hexahedral :
            switch (face)
              {
              case 1 :
                *faceidit = cfxFACEID(elem, 6);
                break;
              case 2 :
                *faceidit = cfxFACEID(elem, 4);
                break;
              case 3 :
                *faceidit = cfxFACEID(elem, 1);
                break;
              case 4 :
                *faceidit = cfxFACEID(elem, 2);
                break;
              case 5 :
                *faceidit = cfxFACEID(elem, 3);
                break;
              case 6 :
                *faceidit = cfxFACEID(elem, 5);
                break;
              default :
                GISMO_ERROR("Invalid face number!");
              }
            break;

          default:
            GISMO_ERROR("Invalid element type!");
          }
      }
  }

  /// \brief Checks if the given element is valid
  static bool isValidElement(const gsCFX5Element& element,
                             const gsCFX5Domain* domain)
  {
    switch (element.getElementType())
      {
      case gsCFX5ElementType::Tetrahedral :
        // Do nothing
        return true;

      case gsCFX5ElementType::Pyramid :
        // Do nothing
        return true;

      case gsCFX5ElementType::Prism :
        // Do nothing
        return true;

      case gsCFX5ElementType::Hexahedral :

        // Check for positive volume using expressions
        // (2.2.29)-(2.2.31) from the lecture notes "NUMERIEKE
        // STROMINGSLEER III" by Prof.dr.ir. P. Wesseling, Delft
        // University of Technology

        double x1[3],x2[3],x3[3],x4[3],x5[3],x6[3],x7[3],x8[3];
        double b1[3],b2[3],b3[3],b12[3],b13[3],b23[3];
        double b2xb3[3],b12xb13[3],b2xb23[3],b23xb3[3];

        x1[0] = domain->getNode(element.getNodeID(0)-1).getX();
        x1[1] = domain->getNode(element.getNodeID(0)-1).getY();
        x1[2] = domain->getNode(element.getNodeID(0)-1).getZ();

        x2[0] = domain->getNode(element.getNodeID(4)-1).getX();
        x2[1] = domain->getNode(element.getNodeID(4)-1).getY();
        x2[2] = domain->getNode(element.getNodeID(4)-1).getZ();

        x3[0] = domain->getNode(element.getNodeID(5)-1).getX();
        x3[1] = domain->getNode(element.getNodeID(5)-1).getY();
        x3[2] = domain->getNode(element.getNodeID(5)-1).getZ();

        x4[0] = domain->getNode(element.getNodeID(1)-1).getX();
        x4[1] = domain->getNode(element.getNodeID(1)-1).getY();
        x4[2] = domain->getNode(element.getNodeID(1)-1).getZ();

        x5[0] = domain->getNode(element.getNodeID(2)-1).getX();
        x5[1] = domain->getNode(element.getNodeID(2)-1).getY();
        x5[2] = domain->getNode(element.getNodeID(2)-1).getZ();

        x6[0] = domain->getNode(element.getNodeID(6)-1).getX();
        x6[1] = domain->getNode(element.getNodeID(6)-1).getY();
        x6[2] = domain->getNode(element.getNodeID(6)-1).getZ();

        x7[0] = domain->getNode(element.getNodeID(7)-1).getX();
        x7[1] = domain->getNode(element.getNodeID(7)-1).getY();
        x7[2] = domain->getNode(element.getNodeID(7)-1).getZ();

        x8[0] = domain->getNode(element.getNodeID(3)-1).getX();
        x8[1] = domain->getNode(element.getNodeID(3)-1).getY();
        x8[2] = domain->getNode(element.getNodeID(3)-1).getZ();

        // Compute formulae (2.2.29)-(2.2.30)
        for (int i=0; i<3; i++)
          {
            b1[i]   = 0.125 * ( x3[i]+x4[i]+x7[i]+x8[i] - x1[i]-x2[i]-x5[i]-x6[i] );
            b2[i]   = 0.125 * ( x2[i]+x3[i]+x6[i]+x7[i] - x1[i]-x4[i]-x5[i]-x8[i] );
            b3[i]   = 0.125 * ( x5[i]+x6[i]+x7[i]+x8[i] - x1[i]-x2[i]-x3[i]-x4[i] );
            b12[i]  = 0.125 * ( x1[i]+x3[i]+x5[i]+x7[i] - x2[i]-x4[i]-x6[i]-x8[i] );
            b13[i]  = 0.125 * ( x1[i]+x2[i]+x7[i]+x8[i] - x3[i]-x4[i]-x5[i]-x6[i] );
            b23[i]  = 0.125 * ( x1[i]+x4[i]+x6[i]+x7[i] - x2[i]-x3[i]-x5[i]-x8[i] );
          }

        // Compute cross-products used in formula (2.2.31)
        b2xb3[0]   = b2[1]*b3[2]   - b2[2]*b3[1];
        b2xb3[1]   = b2[2]*b3[0]   - b2[0]*b3[2];
        b2xb3[2]   = b2[0]*b3[1]   - b2[1]*b3[0];

        b12xb13[0] = b12[1]*b13[2] - b12[2]*b13[1];
        b12xb13[1] = b12[2]*b13[0] - b12[0]*b13[2];
        b12xb13[2] = b12[0]*b13[1] - b12[1]*b13[0];

        b2xb23[0]  = b2[1]*b23[2]  - b2[2]*b23[1];
        b2xb23[1]  = b2[2]*b23[0]  - b2[0]*b23[2];
        b2xb23[2]  = b2[0]*b23[1]  - b2[1]*b23[0];

        b23xb3[0]  = b23[1]*b3[2]  - b23[2]*b3[1];
        b23xb3[1]  = b23[2]*b3[0]  - b23[0]*b3[2];
        b23xb3[2]  = b23[0]*b3[1]  - b23[1]*b3[0];

        double volume;
        volume = 8.0 * (b1[0]*b2xb3[0]   + b1[1]*b2xb3[1]   + b1[2]*b2xb3[2])
          +  8.0/3.0 * (b1[0]*b12xb13[0] + b1[1]*b12xb13[1] + b1[2]*b12xb13[2])
          +  8.0/3.0 * (b12[0]*b2xb23[0] + b12[1]*b2xb23[1] + b12[2]*b2xb23[2])
          +  8.0/3.0 * (b13[0]*b23xb3[0] + b13[1]*b23xb3[1] + b13[2]*b23xb3[2]);

        return (volume < 0);

      default:
        GISMO_ERROR("Invalid element type!");
      }
  }

} // namespace gismo
