SUBROUTINE ONXHalt(String)
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  CHARACTER (LEN=*) :: String
  CALL MondoHalt(EXIT_ERROR,String)
END SUBROUTINE ONXHalt

